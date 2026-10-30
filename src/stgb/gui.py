"""Dark desktop GUI for STGB-Neutron-Stars.

The GUI intentionally delegates all physics calculations to stgb.cli.main() so the
CLI and Windows executable use the same scientific implementation.
"""
from __future__ import annotations

import contextlib
import io
import multiprocessing
import queue
import re
import sys
import threading
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.image as mpimg

from . import cli
from .eos import NAMES

APP_VERSION = "1.3.1"
FIGURE_NAMES = {
    1: "Scalarization mechanism",
    2: "Mass-radius relation",
    3: "Scalar charge",
    4: "Coupling scan",
    5: "Six-EOS comparison",
    6: "Internal structure",
    7: "Linear scalarization threshold",
}


def resource_path(relative: str) -> Path:
    """Resolve repository resources both from source and a PyInstaller bundle."""
    root = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parents[2]))
    return root / relative


class QueueWriter(io.TextIOBase):
    def __init__(self, events: queue.Queue):
        self.events = events
        self._buffer = ""

    def write(self, text: str) -> int:
        if not text:
            return 0
        self._buffer += text
        while "\n" in self._buffer:
            line, self._buffer = self._buffer.split("\n", 1)
            self.events.put(("log", line))
        return len(text)

    def flush(self) -> None:
        if self._buffer:
            self.events.put(("log", self._buffer))
            self._buffer = ""


class STGBApp(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title(f"STGB neutron star — v{APP_VERSION}")
        self.geometry("1220x780")
        self.minsize(1050, 680)
        self.configure(bg="#07090d")

        self.events: queue.Queue = queue.Queue()
        self.running = False
        self.progress_done = 0
        self.progress_total = 1
        self.custom_eos: Path | None = None
        self.last_run: Path | None = None
        self.chart_paths: list[Path] = []

        self._configure_style()
        self._build_ui()
        self.after(100, self._poll_events)

    def _configure_style(self) -> None:
        style = ttk.Style(self)
        try:
            style.theme_use("clam")
        except tk.TclError:
            pass
        bg, panel, fg, muted, accent = "#07090d", "#11151c", "#eef2f7", "#9aa4b2", "#5aa9ff"
        style.configure(".", background=bg, foreground=fg, fieldbackground="#0c1016")
        style.configure("Panel.TFrame", background=panel)
        style.configure("TLabel", background=bg, foreground=fg)
        style.configure("Panel.TLabel", background=panel, foreground=fg)
        style.configure("Muted.TLabel", background=panel, foreground=muted)
        style.configure("Title.TLabel", background=bg, foreground=fg, font=("Segoe UI Semibold", 19))
        style.configure("Section.TLabel", background=panel, foreground=fg, font=("Segoe UI Semibold", 11))
        style.configure("TButton", background="#1a2230", foreground=fg, padding=(10, 7), borderwidth=0)
        style.map("TButton", background=[("active", "#253147"), ("disabled", "#10141b")])
        style.configure("Accent.TButton", background=accent, foreground="#07111d", font=("Segoe UI Semibold", 10))
        style.map("Accent.TButton", background=[("active", "#83beff")])
        style.configure("TCheckbutton", background=panel, foreground=fg)
        style.map("TCheckbutton", background=[("active", panel)])
        style.configure("TCombobox", fieldbackground="#0c1016", background="#0c1016", foreground=fg)
        style.configure("TEntry", fieldbackground="#0c1016", foreground=fg)
        style.configure("Horizontal.TProgressbar", troughcolor="#11151c", background=accent, bordercolor="#11151c", lightcolor=accent, darkcolor=accent)
        style.configure("Treeview", background="#0c1016", fieldbackground="#0c1016", foreground=fg, rowheight=28, borderwidth=0)
        style.map("Treeview", background=[("selected", "#1e4f7a")])

    def _build_ui(self) -> None:
        header = ttk.Frame(self, padding=(22, 18, 22, 10))
        header.pack(fill="x")
        ttk.Label(header, text="STGB neutron star", style="Title.TLabel").pack(side="left")
        ttk.Label(header, text=f"Desktop research interface  •  v{APP_VERSION}", foreground="#9aa4b2").pack(side="left", padx=14, pady=(6, 0))

        main = ttk.Frame(self, padding=(18, 8, 18, 18))
        main.pack(fill="both", expand=True)
        main.columnconfigure(1, weight=1)
        main.rowconfigure(0, weight=1)

        controls = ttk.Frame(main, style="Panel.TFrame", padding=18)
        controls.grid(row=0, column=0, sticky="nsw", padx=(0, 12))
        controls.configure(width=340)
        controls.grid_propagate(False)

        review = ttk.Frame(main, style="Panel.TFrame", padding=12)
        review.grid(row=0, column=1, sticky="nsew")
        review.columnconfigure(0, weight=1)
        review.rowconfigure(1, weight=1)

        self._build_controls(controls)
        self._build_review(review)

    def _build_controls(self, parent: ttk.Frame) -> None:
        ttk.Label(parent, text="Run configuration", style="Section.TLabel").pack(anchor="w")
        ttk.Label(parent, text="The GUI uses the same solver and figure code as the CLI.", style="Muted.TLabel", wraplength=290).pack(anchor="w", pady=(3, 13))

        self.eos_var = tk.StringVar(value="SLy4")
        self.eos_combo = ttk.Combobox(parent, textvariable=self.eos_var, values=NAMES, state="readonly", width=25)
        self._field(parent, "Equation of state", self.eos_combo)
        eos_row = ttk.Frame(parent, style="Panel.TFrame")
        eos_row.pack(fill="x", pady=(0, 10))
        ttk.Button(eos_row, text="Load custom EOS…", command=self._choose_eos).pack(side="left")
        self.custom_eos_label = ttk.Label(eos_row, text="Built-in", style="Muted.TLabel")
        self.custom_eos_label.pack(side="left", padx=8)

        self.alpha_var = tk.StringVar(value="30")
        self.beta_var = tk.StringVar(value="-10")
        self.rho_var = tk.StringVar(value="8.1e14")
        self.quality_var = tk.StringVar(value="quick")
        self._field(parent, "Gauss–Bonnet α [km²]", ttk.Entry(parent, textvariable=self.alpha_var))
        self._field(parent, "Ricci βV", ttk.Entry(parent, textvariable=self.beta_var))
        self._field(parent, "Central ε/c² [g cm⁻³]", ttk.Entry(parent, textvariable=self.rho_var))
        self._field(parent, "Quality", ttk.Combobox(parent, textvariable=self.quality_var, values=("quick", "thesis"), state="readonly"))

        ttk.Label(parent, text="Figures", style="Section.TLabel").pack(anchor="w", pady=(9, 4))
        figures = ttk.Frame(parent, style="Panel.TFrame")
        figures.pack(fill="x")
        self.figure_vars: dict[int, tk.BooleanVar] = {}
        for n, name in FIGURE_NAMES.items():
            var = tk.BooleanVar(value=True)
            self.figure_vars[n] = var
            ttk.Checkbutton(figures, text=f"{n}. {name}", variable=var).pack(anchor="w", pady=1)

        self.output_var = tk.StringVar(value=str(Path.home() / "Documents" / "STGB_Results"))
        ttk.Label(parent, text="Save to", style="Section.TLabel").pack(anchor="w", pady=(13, 5))
        outrow = ttk.Frame(parent, style="Panel.TFrame")
        outrow.pack(fill="x")
        ttk.Entry(outrow, textvariable=self.output_var, width=28).pack(side="left", fill="x", expand=True)
        ttk.Button(outrow, text="…", width=3, command=self._choose_output).pack(side="left", padx=(5, 0))

        self.run_button = ttk.Button(parent, text="Run calculation", style="Accent.TButton", command=self._start_run)
        self.run_button.pack(fill="x", pady=(16, 8))
        self.progress = ttk.Progressbar(parent, mode="determinate", maximum=100)
        self.progress.pack(fill="x")
        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(parent, textvariable=self.status_var, style="Muted.TLabel", wraplength=290).pack(anchor="w", pady=(5, 0))

    def _field(self, parent: ttk.Frame, label: str, widget: tk.Widget) -> None:
        ttk.Label(parent, text=label, style="Muted.TLabel").pack(anchor="w", pady=(6, 3))
        widget.pack(fill="x", pady=(0, 3))

    def _build_review(self, parent: ttk.Frame) -> None:
        top = ttk.Frame(parent, style="Panel.TFrame")
        top.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        ttk.Label(top, text="Chart review", style="Section.TLabel").pack(side="left")
        ttk.Button(top, text="Open results folder", command=self._open_results).pack(side="right")

        body = ttk.Panedwindow(parent, orient="horizontal")
        body.grid(row=1, column=0, sticky="nsew")
        left = ttk.Frame(body, style="Panel.TFrame", padding=(0, 0, 8, 0))
        right = ttk.Frame(body, style="Panel.TFrame")
        body.add(left, weight=1)
        body.add(right, weight=5)

        self.chart_tree = ttk.Treeview(left, show="tree", selectmode="browse", height=16)
        self.chart_tree.pack(fill="both", expand=True)
        self.chart_tree.bind("<<TreeviewSelect>>", self._show_selected_chart)

        self.figure = Figure(figsize=(8, 5), dpi=100, facecolor="#11151c")
        self.ax = self.figure.add_subplot(111)
        self.ax.set_facecolor("#11151c")
        self.ax.text(.5, .5, "Run a calculation to review saved figures here.", ha="center", va="center", color="#9aa4b2", transform=self.ax.transAxes)
        self.ax.set_axis_off()
        self.canvas = FigureCanvasTkAgg(self.figure, master=right)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        toolbar = NavigationToolbar2Tk(self.canvas, right, pack_toolbar=False)
        toolbar.configure(background="#11151c")
        toolbar.pack(fill="x")

        ttk.Label(parent, text="Computation log", style="Section.TLabel").grid(row=2, column=0, sticky="w", pady=(10, 5))
        self.log = tk.Text(parent, height=8, bg="#080b10", fg="#c9d2df", insertbackground="white", relief="flat", padx=9, pady=7, font=("Cascadia Mono", 9))
        self.log.grid(row=3, column=0, sticky="ew")
        self.log.configure(state="disabled")

    def _choose_eos(self) -> None:
        path = filedialog.askopenfilename(title="Select standardized EOS CSV", filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
        if not path:
            return
        self.custom_eos = Path(path)
        self.custom_eos_label.configure(text=self.custom_eos.name)
        self.eos_combo.configure(state="disabled")
        self.figure_vars[5].set(False)
        self.status_var.set("Custom EOS loaded. Figure 5 is disabled because it is the built-in six-EOS comparison.")

    def _choose_output(self) -> None:
        path = filedialog.askdirectory(title="Choose results folder")
        if path:
            self.output_var.set(path)

    def _selected_figures(self) -> list[int]:
        return [n for n, var in self.figure_vars.items() if var.get()]

    def _estimate_work(self, figures: list[int]) -> int:
        quality = self.quality_var.get()
        points = 5 if quality == "quick" else 17
        seq_defaults = {-100.0, 30.0, 200.0} if quality == "quick" else {-200.0, 0.0, 30.0, 100.0, 200.0}
        try:
            seq_defaults.add(float(self.alpha_var.get()))
        except ValueError:
            pass
        total = 0
        if set(figures) & {1, 6}:
            total += 1
        if set(figures) & {2, 3}:
            total += len(seq_defaults) * points
        if 5 in figures:
            total += 6 * points
        if 4 in figures:
            total += 4 if quality == "quick" else 21
        if 7 in figures:
            total += 4 if quality == "quick" else 19
        return max(1, total)

    def _start_run(self) -> None:
        if self.running:
            return
        figures = self._selected_figures()
        if not figures:
            messagebox.showerror("No figures selected", "Select at least one figure.")
            return
        if self.custom_eos and 5 in figures:
            messagebox.showerror("Custom EOS", "Figure 5 is reserved for the six built-in EOS tables. Deselect Figure 5 for a custom EOS run.")
            return
        try:
            float(self.alpha_var.get()); float(self.beta_var.get()); float(self.rho_var.get())
        except ValueError:
            messagebox.showerror("Invalid parameter", "Alpha, beta, and central density must be valid numbers.")
            return

        output = Path(self.output_var.get()).expanduser()
        output.mkdir(parents=True, exist_ok=True)
        before = {p.resolve() for p in output.glob("STGB_*") if p.is_dir()}
        self.progress_done = 0
        self.progress_total = self._estimate_work(figures)
        self.progress["value"] = 1
        self.running = True
        self.run_button.configure(state="disabled")
        self.status_var.set("Calculating…")
        self._clear_log()

        args = [
            "--eos", self.eos_var.get(),
            "--rho-c", self.rho_var.get(),
            "--alpha", self.alpha_var.get(),
            "--beta", self.beta_var.get(),
            "--quality", self.quality_var.get(),
            "--figures", ",".join(map(str, figures)),
            "--data-dir", str(resource_path("data")),
            "--output", str(output),
            "--workers", "1",
            "--no-show", "--save",
        ]
        if self.custom_eos:
            args += ["--eos-file", str(self.custom_eos)]

        threading.Thread(target=self._worker, args=(args, output, before), daemon=True).start()

    def _worker(self, args: list[str], output: Path, before: set[Path]) -> None:
        writer = QueueWriter(self.events)
        try:
            with contextlib.redirect_stdout(writer), contextlib.redirect_stderr(writer):
                rc = cli.main(args)
            writer.flush()
            after = [p for p in output.glob("STGB_*") if p.is_dir() and p.resolve() not in before]
            run_folder = max(after, key=lambda p: p.stat().st_mtime) if after else None
            self.events.put(("done", rc, run_folder))
        except Exception as exc:
            writer.flush()
            self.events.put(("error", repr(exc)))

    def _poll_events(self) -> None:
        try:
            while True:
                event = self.events.get_nowait()
                kind = event[0]
                if kind == "log":
                    self._append_log(event[1])
                    self._update_progress_from_line(event[1])
                elif kind == "done":
                    _, rc, folder = event
                    self.running = False
                    self.run_button.configure(state="normal")
                    self.progress["value"] = 100
                    self.last_run = folder
                    if rc == 0 and folder:
                        self.status_var.set(f"Complete — {folder.name}")
                        self._load_charts(folder)
                    else:
                        self.status_var.set("Calculation finished, but the results folder could not be identified.")
                elif kind == "error":
                    self.running = False
                    self.run_button.configure(state="normal")
                    self.status_var.set("Calculation failed")
                    messagebox.showerror("STGB calculation failed", event[1])
        except queue.Empty:
            pass
        self.after(100, self._poll_events)

    def _update_progress_from_line(self, line: str) -> None:
        counted = False
        if re.search(r"density\s+\d+/\d+", line):
            counted = True
        elif "scanning alpha=" in line or "onset alpha=" in line:
            counted = True
        elif "Figure 1/6: central configuration" in line:
            counted = True
        if counted:
            self.progress_done += 1
            pct = min(96.0, 100.0 * self.progress_done / self.progress_total)
            self.progress["value"] = pct
            self.status_var.set(f"Calculating… {pct:.0f}%")

    def _append_log(self, line: str) -> None:
        self.log.configure(state="normal")
        self.log.insert("end", line + "\n")
        self.log.see("end")
        self.log.configure(state="disabled")

    def _clear_log(self) -> None:
        self.log.configure(state="normal")
        self.log.delete("1.0", "end")
        self.log.configure(state="disabled")

    def _load_charts(self, run_folder: Path) -> None:
        self.chart_paths = sorted((run_folder / "figures").glob("*.png"))
        for item in self.chart_tree.get_children():
            self.chart_tree.delete(item)
        for idx, path in enumerate(self.chart_paths):
            self.chart_tree.insert("", "end", iid=str(idx), text=path.stem.replace("_", " "))
        if self.chart_paths:
            self.chart_tree.selection_set("0")
            self._display_image(self.chart_paths[0])

    def _show_selected_chart(self, _event=None) -> None:
        selected = self.chart_tree.selection()
        if not selected:
            return
        idx = int(selected[0])
        if 0 <= idx < len(self.chart_paths):
            self._display_image(self.chart_paths[idx])

    def _display_image(self, path: Path) -> None:
        image = mpimg.imread(path)
        self.ax.clear()
        self.ax.imshow(image)
        self.ax.set_axis_off()
        self.figure.patch.set_facecolor("#11151c")
        self.canvas.draw_idle()

    def _open_results(self) -> None:
        path = self.last_run or Path(self.output_var.get()).expanduser()
        if not path.exists():
            messagebox.showinfo("Results", "No results folder exists yet.")
            return
        try:
            import os
            if sys.platform.startswith("win"):
                os.startfile(path)  # type: ignore[attr-defined]
            elif sys.platform == "darwin":
                import subprocess; subprocess.Popen(["open", str(path)])
            else:
                import subprocess; subprocess.Popen(["xdg-open", str(path)])
        except Exception as exc:
            messagebox.showerror("Open folder", str(exc))


def main() -> int:
    multiprocessing.freeze_support()
    app = STGBApp()
    app.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

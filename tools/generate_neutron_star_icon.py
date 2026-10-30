"""Generate a simple neutron-star .ico used by the Windows executable build."""
from pathlib import Path
from PIL import Image, ImageDraw, ImageFilter

SIZE = 256
out = Path(__file__).resolve().parents[1] / "assets" / "neutron_star.ico"
out.parent.mkdir(parents=True, exist_ok=True)

img = Image.new("RGBA", (SIZE, SIZE), (4, 6, 10, 255))
glow = Image.new("RGBA", img.size, (0, 0, 0, 0))
g = ImageDraw.Draw(glow)
for radius, alpha in ((92, 20), (70, 34), (52, 62)):
    box = (128-radius, 128-radius, 128+radius, 128+radius)
    g.ellipse(box, fill=(80, 165, 255, alpha))
glow = glow.filter(ImageFilter.GaussianBlur(18))
img = Image.alpha_composite(img, glow)

d = ImageDraw.Draw(img)
d.ellipse((67, 67, 189, 189), fill=(18, 50, 86, 255), outline=(135, 207, 255, 255), width=4)
d.arc((39, 92, 217, 164), start=194, end=346, fill=(210, 231, 255, 255), width=8)
d.arc((46, 79, 210, 178), start=14, end=166, fill=(89, 166, 226, 255), width=4)
d.ellipse((106, 95, 138, 127), fill=(225, 246, 255, 255))
d.ellipse((116, 104, 130, 118), fill=(255, 255, 255, 255))

img.save(out, format="ICO", sizes=[(16,16), (24,24), (32,32), (48,48), (64,64), (128,128), (256,256)])
print(out)

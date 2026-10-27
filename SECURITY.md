# Security Policy

## Supported versions

This is a small research project rather than a security-maintained production service. Security fixes, when possible, will generally target the latest code on the `main` branch. Older revisions and experimental branches are not guaranteed to receive backported fixes or security updates.

## Reporting a vulnerability

**Please do not open a public GitHub issue with details of an exploitable security vulnerability.**

Preferred method: use GitHub's **Report a vulnerability** feature on this repository's **Security** tab, if private vulnerability reporting has been enabled by the maintainer.

If private vulnerability reporting is unavailable, contact the maintainer privately through a contact method listed on the maintainer's GitHub profile. Do not disclose exploit details publicly before the maintainer has had a reasonable opportunity to assess them.

Include, where possible:

- A concise description of the vulnerability and potential impact.
- A minimal reproduction or proof of concept that does not expose private data.
- The affected version, commit, operating system, and relevant dependencies.
- A suggested mitigation or fix, if known.

Reports will be reviewed on a best-effort basis. No fixed response or remediation timeframe is guaranteed.

## Scope and safe use

Potential security issues include unsafe handling of external files, inadvertent exposure of private information, dangerous dependency behavior, or command execution vulnerabilities. Scientific inaccuracies, convergence failures, and incorrect physical predictions should normally be reported through the public numerical/physics issue template, unless they also create a security impact.

Only run code and load EOS datasets from sources you trust. Review third-party dependencies and scripts before executing them. The software is provided **as is**, without a guarantee of fitness for safety-critical, clinical, or operational decision-making.

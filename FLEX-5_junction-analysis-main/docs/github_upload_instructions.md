# GitHub upload instructions

## 1. Create a new repository on GitHub

Create a new repository named `FLEX-5-junction-analysis`. Do not initialize it with a README, license, or `.gitignore`.

## 2. Initialize Git locally

```bash
cd "FLEX-5'junction_analysis"
git init
git add .
git commit -m "Initial commit: FLEX 5' junction analysis pipeline"
```

## 3. Connect and push

HTTPS:

```bash
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/FLEX-5-junction-analysis.git
git push -u origin main
```

SSH:

```bash
git branch -M main
git remote add origin git@github.com:YOUR_USERNAME/FLEX-5-junction-analysis.git
git push -u origin main
```

## 4. Verify after upload

```bash
git clone https://github.com/YOUR_USERNAME/FLEX-5-junction-analysis.git
cd FLEX-5-junction-analysis
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
export FLEX_FASTQ_DIR=/path/to/your/fastq
python analyze_flex.py
```

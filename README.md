# HydroMorph — Web App

**Browser-based hydrocephalus morphometrics pipeline.**
Evans Index · Callosal Angle · Ventricle Volume · NPH Scoring

100% client-side. No data ever leaves your device.

## Deploy (2 steps)

1. **Create a new GitHub repository** and push these files to `main`
2. **Enable GitHub Pages**: Settings → Pages → Source: **GitHub Actions**

That's it. The workflow deploys automatically on every push.

Your app will be live at: `https://<username>.github.io/<repo-name>/`

## What happens on push

```
push to main
  → GitHub Actions triggers
  → Uploads static files (HTML/CSS/JS)
  → Deploys to GitHub Pages
  → Live URL in ~60 seconds
```

## Files

| File | Size | Purpose |
|------|------|---------|
| `index.html` | 13 KB | 3-screen app shell (Upload → Processing → Results) |
| `style.css` | 22 KB | GitHub-dark responsive theme |
| `app.js` | 42 KB | Full pipeline: segmentation, morphometrics, rendering |
| `nifti-reader.js` | 7 KB | NIfTI-1 parser with gzip/endianness |
| `sample-data.json` | 430 KB | Bundled 64×64 CT for "Try Sample" button |

## Manual trigger

You can also trigger a deploy manually:
Actions tab → "Deploy to GitHub Pages" → "Run workflow"

## Author

**Matheus Machado Rech**

Research use only — not for clinical diagnosis.

## License

Data reference: CADS BrainCT-1mm (CC BY 4.0)

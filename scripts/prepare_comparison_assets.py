"""Crop existing vector PDFs for side-by-side review; never recompute physics."""
from pathlib import Path
import hashlib
import json
import pymupdf

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'slides/assets/comparison'
OUT.mkdir(parents=True, exist_ok=True)
records = []

def crop(source, name, rect=None):
    path = ROOT / source
    with pymupdf.open(path) as doc:
        clip = pymupdf.Rect(rect) if rect else doc[0].rect
        result = pymupdf.open()
        page = result.new_page(width=clip.width, height=clip.height)
        page.show_pdf_page(page.rect, doc, 0, clip=clip)
        result.save(OUT / (name + '.pdf'), garbage=4, deflate=True)
        result.close()
    records.append(dict(asset=name+'.pdf', source=source,
                        source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                        page=1, clip_pdf_points=list(clip)))

paper = 'external_data/paper/source/plots/'
fig = 'outputs/figures/'
crop(paper+'IMF_plot.pdf', 'fig1_original')
crop(paper+'tquench_exs.pdf', 'fig2_original')
for z, left, right in [
    (3, (0,0,660,531), (0,20,316,244.3572)),
    (6, (662,0,1243.0455,531), (317,20,583,244.3572)),
    (9, (291,535,952,1072.1149), (584,20,850.421,244.3572)),
]:
    crop(paper+'mh_max_tmax.pdf', f'fig4_z{z}_original', left)
    crop(fig+'figure4_conventions.pdf', f'fig4_z{z}_computed', right)
for key, left, right in [
    ('eff_halo',(0,0,490,430),(0,19,350,290.8)),
    ('eff_stellar',(495,0,914.3079,430),(350,19,705.321,290.8)),
    ('imf_halo',(0,430,490,893.7311),(0,294,350,562.3791)),
    ('imf_stellar',(495,430,914.3079,893.7311),(350,294,705.321,562.3791)),
]:
    crop(paper+'tquench_IMFs.pdf', f'fig5_{key}_original', left)
    crop(fig+'figure5_efficiency_imf_sn_only.pdf', f'fig5_{key}_computed', right)
for key, left, right in [
    ('f001',(0,0,566,509.9246),(0,20,386,315.4191)),
    ('f01',(572,0,1115.8075,509.9246),(386,20,777.321,315.4191)),
]:
    crop(paper+'tquench_diff_fcovers.pdf', f'fig8_{key}_original', left)
    crop(fig+'figure8_cover_gas_sn_only.pdf', f'fig8_{key}_computed', right)
(OUT / 'manifest.json').write_text(json.dumps(records, indent=2)+'\n')
print(f'{len(records)} vector comparison assets written to {OUT}')

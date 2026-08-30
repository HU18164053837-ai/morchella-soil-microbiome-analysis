from pathlib import Path
from PIL import Image, ImageDraw
root=Path(r"D:\path\to\morchella_microbiome\analysis\manuscript_v1\.render_v3")
for prefix,out in [("ase_fig-","ase_figures_contact.png"),("sum-","summary_contact.png")]:
    pages=sorted(root.glob(prefix+"*.png"),key=lambda p:int(p.stem.split('-')[-1]))
    cellw,cellh=760,1010
    sheet=Image.new('RGB',(cellw*2,cellh*((len(pages)+1)//2)),'#ddd')
    for i,p in enumerate(pages):
        im=Image.open(p).convert('RGB'); im.thumbnail((720,950))
        canvas=Image.new('RGB',(cellw,cellh),'white'); canvas.paste(im,((cellw-im.width)//2,40))
        ImageDraw.Draw(canvas).text((15,10),p.stem,fill='black')
        sheet.paste(canvas,((i%2)*cellw,(i//2)*cellh))
    sheet.save(root/out)


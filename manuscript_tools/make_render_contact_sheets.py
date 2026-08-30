from pathlib import Path
from PIL import Image, ImageOps, ImageDraw

root=Path(r"D:\path\to\morchella_microbiome\analysis\manuscript_v1\.render_ase_v2")
pages=sorted(root.glob("page-*.png"), key=lambda p:int(p.stem.split('-')[1]))
for batch in range(0,len(pages),4):
    ims=[]
    for p in pages[batch:batch+4]:
        im=Image.open(p).convert('RGB')
        im.thumbnail((850,1100))
        canvas=Image.new('RGB',(900,1170),'white')
        canvas.paste(im,((900-im.width)//2,40))
        ImageDraw.Draw(canvas).text((20,10),p.stem,fill='black')
        ims.append(canvas)
    sheet=Image.new('RGB',(1800,2340),'#d8d8d8')
    for i,im in enumerate(ims): sheet.paste(im,((i%2)*900,(i//2)*1170))
    sheet.save(root/f"contact_{batch+1:02d}_{min(batch+4,len(pages)):02d}.png")


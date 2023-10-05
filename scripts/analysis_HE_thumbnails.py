#!/usr/bin/env python

import openslide
from tqdm import tqdm
import glob

for slide_fn in tqdm(glob.glob("data/GLASS_OD/Stainings/H&E Slides/*.ndpi")):
  print(slide_fn)
  with openslide.OpenSlide(slide_fn) as slide:
    scale = 1024 / max(slide.dimensions)
    #print(f"scale: {scale}")
    #print([int(m * scale) for m in slide.dimensions])
    dim_str = "x".join([str(int(m * scale)) for m in slide.dimensions])
    out_fn = "output/figures/HE_thumbnails/" + slide_fn.split("/")[-1].replace(".ndpi","") + "_" + dim_str + "px.png"

    thumb = slide.get_thumbnail([int(m * scale) for m in slide.dimensions])
    #thumb.save(out_fn, format="jpeg")
    thumb.save(out_fn, format="png")



#!/usr/bin/env R

tmp <- read.csv("data/GLASS_OD/Stainings/KI67 - QuPath/measurements/measurements_annotations.tsv", sep="\t")


tmp <- tmp |> 
  dplyr::mutate(Classification = NULL) |> # not really used
  dplyr::mutate(Name = NULL) |> 
  dplyr::rename(qupath_label = Image) |> 
  dplyr::mutate(annotated = grepl("annotated", qupath_label)) |> 
  dplyr::mutate(ki67_staining = trimws(gsub("[ ]?\\[.+\\]","",qupath_label))) |> 
  
  assertr::verify(Object.type == "Annotation") |> 
  dplyr::mutate(Object.type = NULL) |> 

  assertr::verify(Parent == "Root object (Image)") |> 
  dplyr::mutate(Parent = NULL) |> 
  
  assertr::verify(ROI %in% c("Polygon", "Geometry")) |> 
  dplyr::mutate(ROI = NULL) |> 
  
  assertr::verify(!annotated | !is.na(Num.Detections)) |> 
  assertr::verify(!annotated | !is.na(Num.KI67neg)) |> 
  assertr::verify(!annotated | !is.na(Num.KI67pos)) |> 
  assertr::verify(!annotated | !is.na(Num.Other)) |> 

  dplyr::filter(annotated)



tmp.aggregated <- tmp |> 
  dplyr::mutate(`Perimeter.µm` = NULL) |> # not needed
  dplyr::mutate(Object.ID  = NULL) |> 
  dplyr::mutate(Centroid.X.µm = NULL) |> 
  dplyr::mutate(Centroid.Y.µm = NULL) |> 
  dplyr::mutate(qupath_label = NULL) |> 
  dplyr::mutate(annotated = NULL) |> 
  dplyr::group_by(ki67_staining) |> 
  dplyr::summarise(
    Num.Detections = sum(Num.Detections),
    Num.KI67neg = sum(Num.KI67neg),
    Num.KI67pos = sum(Num.KI67pos),
    Num.Other = sum(Num.Other),
    `Area.µm.2` = sum(`Area.µm.2`)
  ) |> 
  
  dplyr::mutate(ki67_pos_per_detected_cells = Num.KI67pos / (Num.KI67pos + Num.KI67neg)) |> 
  dplyr::mutate(ki67_pos_per_area_um2 = Num.KI67pos / `Area.µm.2`)




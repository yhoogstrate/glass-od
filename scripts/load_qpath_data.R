#!/usr/bin/env R

a <- read.csv("data/GLASS_OD/Stainings/KI67 - QuPath/measurements/measurements_annotations.tsv", sep="\t") |> 
  dplyr::mutate(Classification = NULL) |> # not really used
  dplyr::mutate(Name = NULL) |> 
  
  assertr::verify(Object.type == "Annotation" | is.na(Object.type)) |> 
  dplyr::mutate(Object.type = NULL) |> 
  
  assertr::verify(Parent == "Root object (Image)" | is.na(Parent)) |> 
  dplyr::mutate(Parent = NULL) |> 
  
  assertr::verify(!is.na(Num.Other))



sum(!is.na(a$Name))


a |> 
  head()


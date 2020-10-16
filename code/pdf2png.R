# pdf to png

library(pdftools)
w.d <- getwd()

files.all <- list.files(paste0(w.d, '/dag/'))

files <- paste0(paste0(w.d, '/dag/'), grep(".pdf", files.all, value=T))
files.out <- paste0(unlist(strsplit(files, ".pdf")), ".png")
dpi <- rep(400,length(files))
for(i in 1:length(files)){
  pdf_convert(files[i], format = "png", dpi = dpi[i],
              filenames = files.out[i])
}



# model specification diagrams
files.all <- list.files(paste0(w.d, '/model-spec/'))

files <- paste0(paste0(w.d, '/model-spec/'), grep(".pdf", files.all, value=T))
files.out <- paste0(unlist(strsplit(files, ".pdf")), ".png")
dpi <- rep(400,length(files))
for(i in 1:length(files)){
  pdf_convert(files[i], format = "png", dpi = dpi[i],
              filenames = files.out[i])
}



# path diagrams
files.all <- list.files(paste0(w.d, '/path-diagram/'))

files <- paste0(paste0(w.d, '/path-diagram/'), grep(".pdf", files.all, value=T))
files.out <- paste0(unlist(strsplit(files, ".pdf")), ".png")
dpi <- rep(400,length(files))
for(i in 1:length(files)){
  pdf_convert(files[i], format = "png", dpi = dpi[i],
              filenames = files.out[i])
}


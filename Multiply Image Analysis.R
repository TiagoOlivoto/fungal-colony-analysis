library(EBImage)

# =========================================
# 0. Parameters
# =========================================
batch_folder <- "C:/Users/gde267.AD/Downloads/plates_batch/" # folder with plate images
label_image_path <- "C:/Users/gde267.AD/Downloads/label_ref.jpeg" # separate label image for this batch
label_known_area_cm2 <- 4  # known area of the reference label in cm²

# Create overlay folder inside batch folder
overlay_folder <- file.path(batch_folder, "overlay")
if(!dir.exists(overlay_folder)) dir.create(overlay_folder)

# =========================================
# 1. Process the label image once for the batch
# =========================================

# Load & resize label image
label_img <- readImage(label_image_path)
label_img <- resize(label_img, w = 800)  # resize to match plate image width
display(label_img, method="raster")

# Normalize + optional gamma
#   - Usually gamma=1 is fine for neutral label brightness
label_img <- normalize(label_img)
label_img <- label_img ^ 1
display(label_img, method="raster")

# Convert to grayscale
label_gray <- channel(label_img, "gray")
display(label_gray, method="raster")

# Threshold to isolate label
#   - If label is darker than background, use < 
#   - If label is lighter than background, use >
thr <- otsu(label_gray)
mask_lt <- label_gray < thr
mask_gt <- label_gray > thr
label_bin <- if(mean(mask_lt) < mean(mask_gt)) mask_lt else mask_gt
label_bin <- fillHull(label_bin)
display(label_bin, method="raster")

# Morphological cleaning
#   - erode(): removes small specks
#   - dilate(): smooth edges
#   - Adjust brush size and shape depending on label size
label_bin <- erode(label_bin, makeBrush(5, shape="disc"))
label_bin <- dilate(label_bin, makeBrush(5, shape="disc"))
display(label_bin, method="raster")

# Label connected components in label mask
label_label <- bwlabel(label_bin)
label_features <- computeFeatures.shape(label_label)

# Assume largest object is the reference label
largest_label <- which.max(label_features[,"s.area"])
label_area_pix <- label_features[largest_label, "s.area"]
cat("Reference label area in pixels:", label_area_pix, "\n")

# Compute pixel → cm² factor
pixel_to_cm2 <- label_known_area_cm2 / label_area_pix
cat("Pixels → cm² factor:", pixel_to_cm2, "\n")

# =========================================
# 2. Process all plate images
# =========================================
plate_files <- list.files(batch_folder, pattern="\\.jpeg$|\\.jpg$|\\.png$", full.names = TRUE)

# Initialize results table
results <- data.frame(
  image_name = character(),
  colony_area_pix = numeric(),
  label_area_pix = numeric(),
  label_area_cm2 = numeric(),
  colony_area_cm2 = numeric(),
  stringsAsFactors = FALSE
)

for (image_path in plate_files) {
  cat("Processing:", image_path, "\n")
  
  # ---- Load & resize plate image ----
  my_image <- readImage(image_path)
  img_resized <- resize(my_image, w = 800)
  
  # ---- Normalize + gamma correction ----
  img_norm <- normalize(img_resized)
  
  # Gamma controls brightness/contrast balance:
  #   - If colonies are dark on light agar, try gamma < 1 (e.g. 0.7–0.9)
  #   - If colonies are light on dark agar, try gamma > 1 (e.g. 1.2–1.5)
  #   - Adjust gamma per batch if agar or lighting changes
  img_gamma <- img_norm ^ 1.7
  display(img_gamma, method="raster")
  
  # ---- Convert to grayscale ----
  #   - If colonies are pigmented (red/orange), use "red" channel
  #   - If bluish/greenish, try "blue" or "green"
  #   - If agar itself is colored, pick channel with largest contrast
  img_gray <- channel(img_gamma, "gray")
  display(img_gray, method="raster")
  
  # ---- Threshold & clean binary mask ----
  thr <- otsu(img_gray)
  mask_lt <- img_gray < thr   # colonies dark
  mask_gt <- img_gray > thr   # colonies light
  img_bin <- if(mean(mask_lt) < mean(mask_gt)) mask_lt else mask_gt
  img_bin <- fillHull(img_bin)
  
  # Morphological cleaning:
  #   - erode(): removes small specks, shrinks colonies slightly
  #   - dilate(): grows them back, smooths edges
  #   - Adjust brush size/shape depending on colony size
  img_bin <- erode(img_bin, makeBrush(5, shape="Gaussian"))
  img_bin <- dilate(img_bin, makeBrush(7, shape="Gaussian"))
  display(img_bin, method="raster")
  
  # ---- Label connected components ----
  img_label <- bwlabel(img_bin)
  display(colorLabels(img_label), method="raster")  # visual check
  
  # ---- Compute colony features ----
  colony_features <- computeFeatures.shape(img_label)
  largest <- which.max(colony_features[,"s.area"])  # assumes largest object = colony
  colony_area_pix <- colony_features[largest, "s.area"]
  cat("Colony area in pixels:", colony_area_pix, "\n")
  
  # ---- Mask of largest object ----
  largest_mask <- img_label == largest
  display(largest_mask, method="raster")
  
  # ---- Overlay largest object ----
  overlay <- paintObjects(largest_mask, toRGB(img_gray), col="red")
  display(overlay, method="raster")
  
  # ---- Save overlay image ----
  overlay_filename <- file.path(overlay_folder, paste0(tools::file_path_sans_ext(basename(image_path)), "_overlay.png"))
  writeImage(overlay, overlay_filename)
  
  # ---- Compute colony area in cm² ----
  colony_area_cm2 <- colony_area_pix * pixel_to_cm2
  cat("Colony area in cm²:", colony_area_cm2, "\n")
  
  # ---- Save results to table ----
  results <- rbind(results, data.frame(
    image_name = basename(image_path),
    colony_area_pix = colony_area_pix,
    label_area_pix = label_area_pix,
    label_area_cm2 = label_known_area_cm2,
    colony_area_cm2 = colony_area_cm2,
    stringsAsFactors = FALSE
  ))
}

# =========================================
# 3. Write CSV
# =========================================
csv_file <- file.path(batch_folder, "colony_measurements.csv")
write.csv(results, csv_file, row.names = FALSE)
cat("Results saved to:", csv_file, "\n")

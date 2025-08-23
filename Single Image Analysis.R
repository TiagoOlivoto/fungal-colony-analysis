library(EBImage)

# =========================================
# 0. Parameters
# =========================================
image_path <- "C:/Users/gde267.AD/Downloads/danilo2.jpg"
overlay_folder <- "C:/Users/gde267.AD/Downloads/single_image_overlay/"  # folder to save overlay
label_image_path <- "C:/Users/gde267.AD/Downloads/label_ref.jpeg"  # separate label photo
label_known_area_cm2 <- 4  # known area of the reference label in cm²

# Create overlay folder if it doesn't exist
if(!dir.exists(overlay_folder)) dir.create(overlay_folder)

# =========================================
# 1. Load colony image
# =========================================
my_image <- readImage(image_path)
display(my_image, method="raster")

# =========================================
# 2. Resize
# =========================================
img_resized <- resize(my_image, w = 800)
display(img_resized, method="raster")

# =========================================
# 3. Normalize + gamma correction
# =========================================
img_norm <- normalize(img_resized)

# Gamma controls brightness/contrast balance:
#   - If colonies are **dark on light agar**, try gamma < 1 (e.g. 0.7–0.9) → brightens background, enhances colonies.
#   - If colonies are **light on dark agar**, try gamma > 1 (e.g. 1.2–1.5) → darkens background, enhances colonies.
#   - If agar has **strong color (red/purple)**, gamma alone may not help → consider channel extraction below.
img_gamma <- img_norm ^ 0.7 # Try 2+ if colony is light on dark background and plate has reflection
display(img_gamma, method="raster")

# =========================================
# 4. Convert to grayscale
# =========================================
# By default: grayscale ("gray") mixes all channels.
#   - If colonies are **pigmented** (red/orange), use "red" channel.
#   - If colonies are **bluish/greenish**, try "blue" or "green".
#   - If agar itself is colored, choose the channel where **colony–agar contrast is largest**.
img_gray <- channel(img_gamma, "gray")
display(img_gray, method="raster")

# =========================================
# 5. Threshold & clean binary mask
# =========================================
# Otsu threshold adapts automatically, but:
#   - If colonies are faint, you may need `thr <- otsu(img_gray)*0.9` (more sensitive).
#   - If background noise is detected, try `thr <- otsu(img_gray)*1.1` (more conservative).
# thr <- otsu(img_gray)
# 
# img_bin <- img_gray < thr # If background is light and colonies are darker than clear media use <
# img_bin <- fillHull(img_bin)

# This code automatically chooses which <> is best depending on colony color
thr <- otsu(img_gray)

mask_lt <- img_gray < thr   # colonies dark
mask_gt <- img_gray > thr   # colonies light

img_bin <- if(mean(mask_lt) < mean(mask_gt)) mask_lt else mask_gt
img_bin <- fillHull(img_bin)

# Morphological cleaning:
#   - erode(): removes small specks, shrinks colonies slightly
#   - dilate(): grows them back to smooth edges
#   - Brush size: 
#       3–5 → light cleaning
#       7–10 → stronger cleaning
#   - Shape: 
#       "disc" = round objects
#       "diamond" = sharp edges
#       "Gaussian" = soft effect
img_bin <- erode(img_bin, makeBrush(5, shape="Gaussian"))
img_bin <- dilate(img_bin, makeBrush(7, shape="Gaussian"))
display(img_bin, method="raster")

# =========================================
# 6. Label connected components
# =========================================
img_label <- bwlabel(img_bin)
display(colorLabels(img_label), method="raster")  # visual check

# =========================================
# 7. Compute colony features
# =========================================
colony_features <- computeFeatures.shape(img_label)
largest <- which.max(colony_features[,"s.area"])  # assumes largest object = colony
colony_area_pix <- colony_features[largest, "s.area"]
cat("Colony area in pixels:", colony_area_pix, "\n")

# =========================================
# 8. Mask of largest object
# =========================================
largest_mask <- img_label == largest
display(largest_mask, method="raster")

# =========================================
# 9. Overlay largest object
# =========================================
overlay <- paintObjects(largest_mask, toRGB(img_gray), col="red")
display(overlay, method="raster")

# =========================================
# LABEL CALIBRATION (separate reference label)
# =========================================

# ==== Load & resize label image ====
label_img <- readImage(label_image_path)
label_img <- resize(label_img, w = 800)
display(label_img, method="raster")

# ==== Normalize + optional gamma ====
label_img <- normalize(label_img)
label_img <- label_img ^ 1  # gamma=1, neutral brightness
display(label_img, method="raster")

# ==== Convert to grayscale ====
label_gray <- channel(label_img, "gray")
display(label_gray, method="raster")

# ==== Threshold to isolate label ====
thr <- otsu(label_gray)
mask_lt <- label_gray < thr   # if label is darker than background
mask_gt <- label_gray > thr   # if label is lighter than background
label_bin <- if(mean(mask_lt) < mean(mask_gt)) mask_lt else mask_gt
label_bin <- fillHull(label_bin)
display(label_bin, method="raster")

# ==== Morphological cleaning ====
label_bin <- erode(label_bin, makeBrush(5, shape="disc"))
label_bin <- dilate(label_bin, makeBrush(5, shape="disc"))
display(label_bin, method="raster")

# ==== Label connected components in label mask ====
label_label <- bwlabel(label_bin)
label_features <- computeFeatures.shape(label_label)

# ==== Assume largest object is the reference label ====
largest_label <- which.max(label_features[,"s.area"])
label_area_pix <- label_features[largest_label, "s.area"]
cat("Reference label area in pixels:", label_area_pix, "\n")

# ==== Compute pixel → cm² factor ====
pixel_to_cm2 <- label_known_area_cm2 / label_area_pix
cat("Pixels → cm² factor:", pixel_to_cm2, "\n")

# =========================================
# 10. Convert colony area to cm² using this factor
# =========================================
colony_area_cm2 <- colony_area_pix * pixel_to_cm2
cat("Colony area in cm²:", colony_area_cm2, "\n")

# =========================================
# 11. Save overlay image
# =========================================
overlay_filename <- file.path(overlay_folder, paste0(tools::file_path_sans_ext(basename(image_path)), "_overlay.png"))
writeImage(overlay, overlay_filename)
cat("Overlay image saved to:", overlay_filename, "\n")
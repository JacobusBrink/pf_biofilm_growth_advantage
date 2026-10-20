// Supplementary Figure S7: segmentation outlines (cyan) on the raw SYTO 9
// image, densest slice of each validation stack. Based on the recorded
// overlay_panel.ijm (F26), applied to all six panels.
// Run in Fiji: Plugins > Macros > Run... and select this file.
// Output: outlines/S7_<id>_<carbon>_<time>.tif (one RGB panel per image)
// Montage (done by hand, 2026-10-14): open the panels in the order F26, F24,
// F15, F18, F22, F07; Image > Stacks > Images to Stack (Scale (largest), no
// bicubic interpolation); Make Montage (3 columns x 2 rows, scale 1, border
// 4, no labels) -> outlines/Figure_S7_montage.tif/.png. Column and row
// labels added in Word.

// folder of this macro (validation/final/), so it runs from any location
base = File.getDirectory(getInfo("macro.filepath"));

ids    = newArray("F26", "F24", "F15", "F18", "F22", "F07");
slices = newArray(5, 5, 6, 6, 8, 6);
names  = newArray("glucose_3h", "glycerol_3h", "succinate_3h",
                  "glucose_7h", "glycerol_7h", "succinate_7h");

setBatchMode(true);
setOption("BlackBackground", true);
setOption("ScaleConversions", true);

for (i = 0; i < ids.length; i++) {
  // raw image: one slice
  open(base + "images/" + ids[i] + ".tif");
  stack = getTitle();
  setSlice(slices[i]);
  run("Duplicate...", "title=raw");
  close(stack);

  // labels: same slice -> cell borders -> binary outlines
  open(base + "_labels/" + ids[i] + ".tif");
  stack = getTitle();
  setSlice(slices[i]);
  run("Duplicate...", "title=labels");
  close(stack);
  selectImage("labels");
  run("Find Edges");
  setThreshold(1, 65535);
  run("Convert to Mask");

  // raw: auto contrast, 8-bit
  selectImage("raw");
  run("Enhance Contrast", "saturated=0.35");
  run("8-bit");

  // grey raw + cyan outlines -> RGB
  run("Merge Channels...", "c4=raw c5=labels create");
  run("Stack to RGB");
  rgb = getTitle();

  // scale and 2 um scale bar (burned in, no text)
  run("Properties...", "channels=1 slices=1 frames=1 unit=micron pixel_width=0.0916 pixel_height=0.0916 voxel_depth=0.0916");
  run("Scale Bar...", "width=2 height=4 color=White background=None location=[Lower Right] horizontal bold hide");

  saveAs("Tiff", base + "outlines/S7_" + ids[i] + "_" + names[i] + ".tif");
  run("Close All");
}

setBatchMode(false);
print("Saved " + ids.length + " panels to " + base + "outlines/");

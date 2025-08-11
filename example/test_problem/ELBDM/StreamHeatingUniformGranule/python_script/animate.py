import imageio
import glob

# animation
xy_images = sorted(glob.glob("Data_000*_Slice_z_Dens.png"))

# skip every 10th frame
skip = 5
xy_images = xy_images[::skip]

with imageio.get_writer("animation_z_slice.gif", mode="I", duration=0.2) as writer:
    for filename in xy_images:
        image = imageio.imread(filename)
        writer.append_data(image)

print(f"Done — used {len(xy_images)} frames (skip={skip})")

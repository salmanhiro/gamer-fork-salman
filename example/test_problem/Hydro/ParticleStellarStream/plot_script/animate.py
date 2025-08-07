import imageio
import glob

# animation
xy_images = sorted(glob.glob("Data_000*.png"))
with imageio.get_writer("animation_xy.gif", mode="I", duration=0.2) as writer:
    for filename in xy_images:
        image = imageio.imread(filename)
        writer.append_data(image)

print("Done")

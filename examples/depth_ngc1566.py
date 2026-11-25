from decals_pipeline_class import DecalsPipeline

if __name__ == "__main__":
    pipe = DecalsPipeline(
        layer="ls-dr9",
        galaxy_centers={"NGC1566": (65.0042, -54.9378)},
    )
    pipe.run(["NGC1566"], size=3000, skip_download=True,
             do_plot=True, do_psf_split=False,
             do_save_split=True, do_depth=True, depth_region="south")

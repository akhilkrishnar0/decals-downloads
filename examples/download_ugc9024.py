from decals_pipeline_class import DecalsPipeline

if __name__ == "__main__":
    pipe = DecalsPipeline(
        layer="ls-dr9",
        galaxy_centers={"UGC9024": (211.669, 22.070)},
    )
    pipe.run(["UGC9024"], size=1500, skip_download=True,
             do_plot=False, do_psf_split=False,
             force_download=False, do_save_split=True)

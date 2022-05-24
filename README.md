# MANTRA
Welcome to the home of MANTRA, the "MountAiN glacier Transient snowline Retrieval Algorithm" -- a Google Earth Engine (GEE) toolset to process and evaluate glacier surface classifications and Transient Snowline Altitudes (TSLA) from Landsat and Digital elevation model data.

**At this point the repository contains only the Evaluate Tool to visualize indidivual results. The routine for parallel processing of comprehensive TSLA datasets will be added upon publication of the paper (which is currently in review).**

<img src="https://github.com/cryotools/mantra/blob/main/supplement/mantra-evaluate-screenshot.png" alt="Screenshot of the MANTRA Evaluate Tool">

Please note that results provided by the Evaluate Tool may be flawed, e.g. by cloud cover. 
It is advised to check the results by comparing the classification with the Landsat image by toggling layers before.

FEATURES
- Classify surface materials based on Landsat band-ratios (snow, ice, debris/rock, clouds)
- Create a catalogue of relevant Landsat scenes from multiple missions (Landsat 4 to 8)
- Merge scenes taken at the same date that are covered by the glacier's outline
- Mask pixels with insufficient illumination
- Visualize the classification
- Calculate and display statistics


# Setup

1. Log in to your [GEE account](https://code.earthengine.google.com/).
2. Head over to the "Assets" tab and add a RGI6 shapefile (e.g. the one provided in the "RGI" directory in this repository).
3. Go back to the "Scripts" tab and create two new files called "evaluate" and "core" or similar.
4. Copy the contents of "evaluate.js" and "core.js" from this repo into the respective files in your GEE environment.
5. Edit the evluate script so that the variables "rgi_path" and "core" point at the correct files in your GEE environment (rgi_path referring to the RGI6 shapefile you just loaded to the Assets). 

That's it, you should now be ready to ride the MANTRA.

# Running the MANTRA Evaluate Tool
1. In the "evaluate" script, edit the variables "rgi_id", "start" and "finish" to meet your study subject:
    - For "rgi_id" enter, well the RGI ID of the glacier. In case you used the RGI shapefile provided with this repo, the initial values should work just fine.
    - "start" and "end" refer to period of time to investigate. Make sure to keep it short (max. ~six months) to avoid GEE stalling.
2. Click "Save" and "Run".
3. When GEE has finsihed loading, the dates for which adequate Landsat scenes are available can be selected (combo box on the top right of the map window). Details of the Landsat scenes can be investigated in the "Console" tab.
4. After selecting a date from the list, GEE will start the processing (may take a while).
5. After processing has finished, the classification results will be shown in the map, statistics in the pane on the right-hand side.


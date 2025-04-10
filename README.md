# Cell Data Dashboard
[![DOI](https://zenodo.org/badge/481438728.svg)](https://doi.org/10.5281/zenodo.15107352)

## Pulling Changes
If a change has been made, you will want to open up your terminal and `cd` (change directory) into the path where your CellDataDashboard is (or do the folder dragging thing on MacOS). Once you are in the folder in your terminal, type
```bash
git pull
```
This will *pull* the latest updated version, which you can then run.
## Setup
### Step 1: Install Python
Head over to https://www.python.org/ and grab the latest version of Python (I wrote this in version `3.8` but `3.10` should be fine). Just click Downloads and download the latest source release. ***Please be sure to tick "Add Python to PATH" during / after the installation***.

To test it out after you have installed it, open up a terminal (Terminal on MacOS, Command Prompt / Powershell / Windows Terminal on Windows), type `python` and hit enter. You'll know it's working if you don't get an error. Exit it by typing `exit()` then pressing enter or pressing `Ctrl-D`. 

### Step 2: Download this repository and install dependencies
Click on the green `Code` button in the upper right and then `Download ZIP`, then unzip it somewhere. After you have done that, you're going to to have to do some hacker(wo)man stuff.

First, you're going to want to get the path of the folder you just unzipped.
- On Windows, this is as simple as just clicking the bar at the top of the Windows Explorer when you're in the folder
- On Mac, you can do some fancy stuff like this: https://themacbeginner.com/full-path-on-terminal-by-dragging-a-folderfile-mac-os-x/, which is strange but convenient

Then, in a terminal, type 
```bash
cd whatever/path/you/got/
```
Once you are in the directory, you will then install the dependencies for the project by running
```bash
pip install -r requirements.txt 
```

Wait for that to finish, and the setup should be completed.

### Step 3: Run the dashboard
After you've installed all the dependencies, while still in your terminal, type:
```bash
streamlit run main_page.py
```

That will start the dashboard (a browser tab should pop open automatically).
### Step 4: Using the dashboard
First, the user must specify the absolute path to the folder containing the data outputted from CellProfiler. After specifying a data path, the user must click the Load button and wait for the data to be processed. Following that, the user can select the slide (QC, Experimental Slide 1-6) and tissue (54 samples of FFPET and FFT) they wish to extract and analyze further from the sidebar. Should the user wish to perform spatial analyses across different sub-localizations, they can swap between pages on the sidebar as well.

### Step 5: Exiting the dashboard
Once you're done with the dashboard, go back to the terminal and press `Ctrl-C` to kill it.

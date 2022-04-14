# Cell Data Dashboard
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
streamlit run dashboard.py
```

That will start the dashboard (a browser tab should pop open automatically).
### Step 4: Profit
Once you're done with the dashboard, go back to the terminal and press `Ctrl-C` to kill it.

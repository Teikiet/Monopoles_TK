"""
Parameters used by other files 
"""
from os import path

# For NUR files, can't split up events/sims
NUR = True

# Directory names
EVTS_PATH = "events"
SIMS_PATH = "sims"
SCRIPTS_PATH = "scripts"
DET_PATH = "det_jsons"
CON_PATH = "config"

# Formatting
EVTS_FORMAT = "evt_{0}.hdf5"
SIMS_FORMAT = "sim_{0}.hdf5"
CON_FORMAT = "con_{0}.yaml"

# Remove + sign from filename
NUM_FORMAT = lambda x: f"{x:.3e}".replace("+", "")

# The power raising 10 to
E_MIN = 18
E_MAX = 18

# Number of simulations, logarithmically equidistant
# between 10^E_MIN and 10^E_MAX
NUM_SIMS = 1

# Number of events per simulation
NUM_EVENTS = 100
NUM_EVENTS_PER_FILE = 1

# Detector's json file
DET_NAME = "ara2"
DET_FILE = path.join(DET_PATH, f"{DET_NAME}_detector.json")
DETECTOR = f"/home/tkiet/nuradiomc/{DET_NAME}.json"

#Configuration:
CONFIG = "/home/tkiet/nuradiomc/config_file.yaml"
PARTICLE_NAME = "monopoles"
# Parent Directory names
EVTS_DIR = f"{EVTS_PATH}_n{NUM_FORMAT(NUM_EVENTS)}_{NUM_SIMS}_{E_MIN}to{E_MAX}"
Finish_DIR = f"finish_{SIMS_PATH}_{DET_NAME}_n{NUM_FORMAT(NUM_EVENTS)}_{NUM_SIMS}_{E_MIN}to{E_MAX}"
SIMS_DIR = f"{SIMS_PATH}_{PARTICLE_NAME}_{DET_NAME}_n{NUM_FORMAT(NUM_EVENTS)}_{NUM_SIMS}_{E_MIN}to{E_MAX}"
SCRIPTS_DIR = f"{SCRIPTS_PATH}_{PARTICLE_NAME}_{DET_NAME}_n{NUM_FORMAT(NUM_EVENTS)}_{NUM_SIMS}_{E_MIN}to{E_MAX}"
CON_DIR = f"{CON_PATH}_{DET_NAME}_n{NUM_FORMAT(NUM_EVENTS)}_{NUM_SIMS}_{E_MIN}to{E_MAX}"


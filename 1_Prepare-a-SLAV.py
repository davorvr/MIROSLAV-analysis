# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: .venv
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Prepare-a-SLAV

# %% [markdown]
# Prepare-a-SLAV utilises the mirofile library to load raw MIROSLAV data into a pandas data frame. Appropriate animal IDs are applied, and the data is downsampled from the default MIROSLAV sampling rate (10 sensor readings/binary values per second) to an arbitrary, user-defined time interval (bin). Prepare-a-SLAV's configuration is performed through the Prepare-a-SLAV TOML configuration file where you can find more information about its parameters.
#
# If you are running Prepare-a-SLAV via Google Colab, Prepare-a-SLAV will autodetect and set up the Colab environment in the following cell, and pull example data and the TOML configuration file from the [MIROSLAV toolkit GitHub repository](https://github.com/davorvr/MIROSLAV-analysis).
#
# If you want to run Prepare-a-SLAV in Google Colab *and* with your own data, you can upload your configuration and data files using the File Browser in the sidebar on the left after running the following cell.

# %%
import importlib.util
try:
    importlib.util.find_spec("google.colab")
except ModuleNotFoundError:
    IN_COLAB = False
else:
    IN_COLAB = True
    import sys
    from IPython.display import clear_output 
    import importlib.metadata
    import packaging.version
    if packaging.version.Version(importlib.metadata.version("pandas")) < packaging.version.Version("2.2"):
        # %pip install -Uq "pandas==2.2"
    # %pip install --ignore-requires-python mirofile
    # %pip install fastparquet
    # !wget https://raw.githubusercontent.com/davorvr/MIROSLAV-analysis/main/1_Prepare-a-SLAV_config.toml
    # !mkdir 0_raw_logs
    # !wget -O 0_raw_logs/mph-pir-rack_M.2022-05-06T19-19-57-669585.gz https://github.com/davorvr/MIROSLAV-analysis/raw/main/0_raw_logs/mph-pir-rack_M.2022-05-06T19-19-57-669585.gz
    # !wget -O 0_raw_logs/mph-pir-rack_M.2022-05-06T19-19-57-669585.gz https://github.com/davorvr/MIROSLAV-analysis/raw/main/0_raw_logs/mph-pir-rack_M.2022-05-06T19-19-57-669585.gz
    # !wget -O 0_raw_logs/mph-pir-rack_M.2022-05-16T01-33-25-478055.gz https://github.com/davorvr/MIROSLAV-analysis/raw/main/0_raw_logs/mph-pir-rack_M.2022-05-16T01-33-25-478055.gz
    # !wget -O 0_raw_logs/mph-pir-rack_M.2022-05-25T14-03-17-240158.gz https://github.com/davorvr/MIROSLAV-analysis/raw/main/0_raw_logs/mph-pir-rack_M.2022-05-25T14-03-17-240158.gz
    # !wget -O 0_raw_logs/mph-pir-rack_R.2022-05-06T19-19-57-669185.gz https://github.com/davorvr/MIROSLAV-analysis/raw/main/0_raw_logs/mph-pir-rack_R.2022-05-06T19-19-57-669185.gz
    # !wget -O 0_raw_logs/mph-pir-rack_R.2022-05-16T01-33-22-935712.gz https://github.com/davorvr/MIROSLAV-analysis/raw/main/0_raw_logs/mph-pir-rack_R.2022-05-16T01-33-22-935712.gz
    # !wget -O 0_raw_logs/mph-pir-rack_R.2022-05-25T07-57-01-575482.gz https://github.com/davorvr/MIROSLAV-analysis/raw/main/0_raw_logs/mph-pir-rack_R.2022-05-25T07-57-01-575482.gz
    clear_output()

# %% [markdown]
# The environment has been set up. If you wish, you can load your own data using the sidebar now.

# %% [markdown]
# ***

# %%
import pandas as pd
if IN_COLAB and sys.hexversion < 0x030b0000:
    OLD_TOML = True
    import toml
else:
    OLD_TOML = False
    import tomllib
import os
from datetime import datetime
from pathlib import Path
from math import ceil
from mirofile import mirofile

# %% [markdown]
# Define helper functions for managing imported column mappings

# %%
def _unpack_mcp_colmap(colmap: dict):
    unpacked_colmap = {}
    # first do PH7..0
    for i in range(7,-1,-1):
        unpacked_colmap.update({"PH"+str(i):"H_"+colmap["PHL"+str(i)]})
    # then do PL0..7
    for i in range(0,8):
        unpacked_colmap.update({"PL"+str(i):"L_"+colmap["PHL"+str(i)]})
    return unpacked_colmap

def unpack_full_colmap(colmap: dict[dict]):
    unpacked_colmap = {}
    input_boards = list(colmap.keys())
    input_boards.remove("top_board")
    input_boards.sort()
    input_boards = ["top_board"] + input_boards
    for mcp_i, colmap_name in enumerate(input_boards):
        mcp_colmap = colmap[colmap_name]
        mcp_colmap_unpacked = _unpack_mcp_colmap(mcp_colmap)
        # new_mcp_colmap = {}
        for k, v in mcp_colmap_unpacked.items():
            #new_mcp_colmap.update({f"MCP{mcp_i+1}_"+k : v})
            unpacked_colmap.update({f"MCP{mcp_i+1}_"+k : v})
        #unpacked_colmap.append(new_mcp_colmap.copy())
    return unpacked_colmap


# %% [markdown]
# Set the current working directory to the location of this script

# %%
wd = Path(os.path.dirname(os.path.realpath('__file__'))).resolve()

# %% [markdown]
# Load the TOML config file and extract all user-defined parameters

# %%
if OLD_TOML:
    with open(wd / "1_Prepare-a-SLAV_config.toml", "r") as cfg_file:
        config = toml.load(cfg_file)
else:
    with open(wd / "1_Prepare-a-SLAV_config.toml", "rb") as cfg_file:
        config = tomllib.load(cfg_file)

try:
    experiment = config["id_variables"]["experiment"]
    set_dtypes = config["processing_params"]["set_dtypes"]
    resample = config["processing_params"]["resample"]
    resample_bin = config["processing_params"]["resample_bin"]
    toml_colnames = config[experiment]
except KeyError as exc:
    raise KeyError("Config file is improperly formatted!") from exc

colmaps = {}
for k, v in toml_colnames.items():
    try:
        do_process = v.pop("process")
    except KeyError as exc:
        raise KeyError("Config file is improperly formatted!") from exc
    if not do_process:
        continue
    try:
        v = unpack_full_colmap(v)
    except Exception as exc:
        raise KeyError("Couldn't process cage mappings from config file!") from exc
    colmaps.update({k: v})

log_path = Path.cwd() / "0_raw_logs"
output_path = Path.cwd() / "1_outputs_prepared"
output_path.mkdir(exist_ok=True)

# %% [markdown]
# The number of rows read at once. Reduce if you run into memory issues.

# %%
chunk_size = 10**6

# %%
ts_column = "ts_recv"
ts_delta = True
ts_index_map = {"ts_recv": 0,
                "ts_sent": 1}
ts_index = ts_index_map[ts_column]
resample_bin_td = pd.to_timedelta(resample_bin)

for device, colmap in colmaps.items():
    print(f"Processing device {device}. ")
    mname = "-".join([experiment, "pir", device])

    mfile = mirofile.open_experiment(experiment, device, compression="gz", path=log_path)
    suffix = ""
    if set_dtypes:
        suffix += "-dtyped"
    if resample:
        suffix += f"-resampled-{resample_bin}".replace(" ", "")
    pqfile = Path(output_path, mname+suffix+".parquet")
    
    if pqfile.exists():
        raise FileExistsError("Database file already exists, not overwriting!")
        #pass

    # Column names for each of the 16 bits outputted by one MCP are stored
    # in order in the mirofile.mcp_colnames list (for more explanation, see
    # mirofile_columns.txt). But, since we usually have more than one MCP,
    # we need to prepend the column names with "MCPn_". This line does this:
    data_columns = [ f"MCP{n_mcp}_"+colname for n_mcp in range(1, mfile.n_mcps+1) for colname in mirofile.mcp_columns ]
    all_columns = mirofile.timestamp_columns+data_columns
    all_column_dtypes = [*["datetime64[ns]"]*2, *["int"]*(len(all_columns)-2)]

    start = datetime.now()
    file_chunk_carryover = None
    #file_chunk = mfile.readlists(size=chunk_size, progress_bar=True)
    sampling_period = None
    n_chunk = 0
    while True:
        n_chunk += 1
        print(f"({device}) Processing chunk {n_chunk}...")
        if isinstance(file_chunk_carryover, pd.DataFrame):
            #file_chunk = [file_chunk_carryover]
            file_chunk = mfile.readlists(size=chunk_size-len(file_chunk_carryover), progress_bar=True)
        else:
            file_chunk = mfile.readlists(size=chunk_size, progress_bar=True)
        if not file_chunk:
            break
        file_chunk = pd.DataFrame.from_records(file_chunk, columns=all_columns)
        if set_dtypes:
            file_chunk = file_chunk.astype(dict(zip(all_columns, all_column_dtypes)))
            if isinstance(file_chunk_carryover, pd.DataFrame) and not file_chunk_carryover.empty:
                file_chunk = pd.concat([file_chunk_carryover, file_chunk], ignore_index=True)
                file_chunk_carryover = None
            if not sampling_period:
                sampling_period = file_chunk[ts_column].diff().mode().iloc[0]
                supp_len = ceil(resample_bin_td / sampling_period)
            if resample:
                file_chunk_supplement = []
                bin_end = file_chunk[ts_column].iloc[-1].ceil(freq=resample_bin_td)
                file_chunk_supplement = mfile.readlists(supp_len, progress_bar=True)
                if file_chunk_supplement:
                    #DEBUG: if pd.to_datetime(last_line[ts_index]) > pd.Timestamp('2022-05-17 05:59:59.000000'):
                        #DEBUG: print("break here")
                    chunk_end = pd.to_datetime(file_chunk_supplement[-1][ts_index])
                    while chunk_end < bin_end:
                        file_chunk_supplement.append(mfile.readlists(supp_len, progress_bar=True))
                        chunk_end = pd.to_datetime(file_chunk_supplement[-1][ts_index])
                        #last_line = mfile.readlist()
                
                    file_chunk_supplement = pd.DataFrame.from_records(file_chunk_supplement, columns=all_columns)
                    file_chunk_supplement = file_chunk_supplement.astype(dict(zip(all_columns, all_column_dtypes)))
                    file_chunk = pd.concat([file_chunk, file_chunk_supplement.loc[file_chunk_supplement[ts_column] < bin_end].copy()], ignore_index=True)
                    file_chunk_carryover = file_chunk_supplement.loc[file_chunk_supplement[ts_column] >= bin_end].copy()
                    file_chunk_supplement = None
                else:
                    file_chunk_carryover = None

        if colmap:
            file_chunk = file_chunk.rename(colmap, axis="columns")
            file_chunk = file_chunk.loc[:,~file_chunk.columns.str.startswith("H_na")]
            file_chunk = file_chunk.loc[:,~file_chunk.columns.str.startswith("L_na")]
            file_chunk = file_chunk.copy()
        if ts_delta:
            if ts_column == "ts_recv":
                secondary_ts = "ts_sent"
            elif ts_column == "ts_sent":
                secondary_ts = "ts_recv"
            else:
                raise ValueError("ts_column must be either 'ts_recv' or 'ts_sent'")
            # delta is always ts_recv - ts_sent since ts_recv is always later, but
            # we keep the other's name so it's clear which column got replaced with a delta.
            file_chunk[secondary_ts+"_delta"] = file_chunk["ts_recv"] - file_chunk["ts_sent"]
            file_chunk = file_chunk.drop(columns=secondary_ts)
        file_chunk = file_chunk.set_index(ts_column)
        if resample:
            if ts_delta:
                file_chunk = file_chunk.resample(resample_bin_td).mean(numeric_only=False)
            else:
                file_chunk = file_chunk.resample(resample_bin_td).mean(numeric_only=True)
            
        if not pqfile.exists():
            file_chunk.to_parquet(pqfile.absolute(), engine="fastparquet")
        else:
            file_chunk.to_parquet(pqfile.absolute(), engine="fastparquet", append=True)
        #DEBUG: last_chunk = file_chunk.copy()
        print(f"({device}) Chunk processed. Last timestamp: {file_chunk.index[-1]}")
    print(f"{pqfile.name}: ", datetime.now()-start)

# %% [markdown]
# ***

# %% [markdown]
# You should be able to obtain the output files from the sidebar now, and proceed to [TidySLAV](https://github.com/davorvr/MIROSLAV-analysis).

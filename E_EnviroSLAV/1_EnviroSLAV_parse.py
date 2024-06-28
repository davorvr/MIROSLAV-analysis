import gzip
import pandas as pd


## Environmental data from a real experiment recorded with a prototype MIROSLAV.
## Load, parse to a dataframe, and save to a .parquet file.
df = pd.read_csv("0_raw_env/0_proto_env.csv.gz",
                 compression="gzip",
                 names=["ts",
                        "presence",
                        "illumination",
                        "temperature",
                        "humidity"],
                 index_col="ts",
                 dtype={"presence": "int64",
                        "illumination": "float64",
                        "temperature": "float64",
                        "humidity": "float64"},
                 parse_dates=True)
df.to_parquet("1_parsed_env/1_proto_env.parquet")

time_bin = "15 minutes"
df = df.resample(pd.to_timedelta(time_bin)).agg({"presence": "max",
                                                 "illumination": "mean",
                                                 "temperature": "mean",
                                                 "humidity": "mean"})
df.to_parquet(f"1_parsed_env/1_proto_env_resampled-{time_bin.replace(" ", "-")}.parquet")


## Test environmental data from a recent MIROSLAV v0.4, sample to illustrate the logging format.
## Load, parse to a dataframe, and save to a .parquet file.
with gzip.open("0_raw_env/miroslav-env-rack_M.2024-06-27T14-04-44-916883.gz", "rt") as f:
    envlog = f.readlines()
for n, line in enumerate(envlog):
    ts_recv, line = line.strip().split(";")
    ts_recv = ts_recv.strip('"')
    if not (line.startswith('"START ') and line.endswith(', END"')):
        raise ValueError(f"Line {n} seems improperly formatted.")
    
    line = line[7:-6].split(",")

    ts_sent = line.pop(0)
    pir, lx, temp, rh = [ x.split(":")[1] for x in line]

    rh = rh.rstrip("%")

    envlog[n] = [ts_recv, ts_sent, pir, lx, temp, rh]

df = pd.DataFrame.from_records(envlog, columns=["ts_recv",
                                                "ts_sent",
                                                "presence",
                                                "illumination",
                                                "temperature",
                                                "humidity"])
df["ts_recv"] = pd.to_datetime(df["ts_recv"])
df["ts_sent"] = pd.to_datetime(df["ts_sent"])
df = df.astype({"presence": "bool",
                "illumination": "float64",
                "temperature": "float64",
                "humidity": "float64"})
df.to_parquet("1_parsed_env/miroslav-env-rack_M.2024-06-27T14-04-44-916883.parquet")
df = df.drop(columns="ts_recv")
df = df.set_index("ts_sent")
df.plot()

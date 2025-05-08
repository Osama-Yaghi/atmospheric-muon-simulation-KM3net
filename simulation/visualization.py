import matplotlib
matplotlib.use("Agg")
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import os

def store_results(results,t0=0,t1=100):
    filtered = [r for r in results if r is not None]
    df = pd.DataFrame(filtered, columns=["x", "y", "z", "vx", "vy", "vz", "E", "num_calls"])

    run_id = int(t1)
    os.makedirs("output", exist_ok=True)
    out_file = f"output/result_{run_id}.csv"
    df.to_csv(out_file, index=False)

    print(f"Simulation complete. {len(df)} muons detected. Output saved to {out_file}.")
    print("Generating plots...")

    fig1,ax1=plt.subplots(2,1,figsize=(20,8))
    sns.histplot(df["E"], bins='auto', kde=False, ax=ax1[0], element="step")
    sns.histplot(df["vz"], bins='auto', kde=False, ax=ax1[1], element="step")
    ax1[0].set_xlabel("Energy (GeV)")
    ax1[0].set_ylabel("Count")
    ax1[0].set_title("Muon Energy Distribution")
    fig1.savefig(f"output/hists1.png")

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(df["x"], df["y"], df["z"],marker='o')
    ax.set_title("Muon Impact Locations")
    plt.savefig(f"output/3d1.png")

    print("Plots saved. Done in {:.2f} seconds.".format(t1 - t0))
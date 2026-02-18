import os
import sys
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from functions.utils import log_message, check_dir

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))


def main():
    log_message("Loading HAR-TF pairs scores...")

    tf_data_path = os.path.join(
        PROJECT_ROOT, "results/har_tf/human/har_tf_pairs_scores.csv"
    )

    if not os.path.exists(tf_data_path):
        log_message(
            f"HAR-TF pairs scores file not found: {tf_data_path}", message_type="error"
        )
        return

    tf_data = pd.read_csv(tf_data_path)

    tfs = tf_data["TF"].unique()
    tfs = sorted(tfs)

    log_message(f"Found {len(tfs)} unique TFs")

    output_dir = check_dir(os.path.join(PROJECT_ROOT, "results/har_tf"))
    output_path = os.path.join(output_dir, "tfs.csv")

    tfs_df = pd.DataFrame(tfs, columns=["TF"])
    tfs_df.to_csv(output_path, index=False, quoting=0)

    log_message(f"Saved unique TFs to: {output_path}", message_type="success")


if __name__ == "__main__":
    main()

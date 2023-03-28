import pandas as pd
import argparse
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parse csvs")
    parser.add_argument(
        '--save_filename', help='Result filename', default="KSK_results", required=False)
    parser.add_argument(
        '--csvs', help='List of CSVs', required=True, default=[], nargs=24)
    args = parser.parse_args()
    print("Toto je args.csvs")
    print(args.csvs)
    result = pd.DataFrame()
    for csv_path in args.csvs:
        csv_path = csv_path.replace("[", "").replace(",", "").replace("]", "")
        print(csv_path)
        if os.stat(csv_path).st_size == 0:
            print("OS stat 0")
            continue
        csv = pd.read_csv(csv_path)
        if csv.empty:
            print("CSV empty")
            continue
        print('Concat ' + csv_path)
        result = pd.concat([result, csv], ignore_index=True)
    result.to_csv(f"{args.save_filename}.csv")

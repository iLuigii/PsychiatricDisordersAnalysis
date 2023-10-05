import os
import json
import sys

# data module is used to retrieve data and process it using the pipeline
import src.data.make_dataset as dm
import src.features.build_features as bf

def main(target):
    if target == "run.py":
        dm.main()
        bf.main()
        print("run the whole project")
    elif target == "test":
        print("start test on making datasets")
        dm.test()
        print("finished running on the test target")

        print("start test on making features\n")
        bf.test()
        print("\nfinished")

    elif target == "check":
        print("start performing checking on the processed data")
        dm.check()
        print("checking process has been performed")
    else:
        print("The command you enter is not correct")

if __name__ == "__main__":
    target = sys.argv[-1]
    main(target)
    
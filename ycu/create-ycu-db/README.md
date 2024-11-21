
# Description
**Create YCU database.**  
:warning: For YCU lab.  
This tiny program is to summarize sample information as a SQLite database.


# Usage
### 1. Copy three Excel files  
Copy the following three Excel files to the `rawdata` directory:
1. Excel file for the old REP
2. Excel file for the current REP being updated
3. Excel file of the mailed list.

### 2. Setup python environment 
Setup the environment using the following commands in a conda-enabled setting.
```Shell
$ conda create -n ycudb python=3.9 -y && conda activate ycudb
$ conda install -y pandas openpyxl sqlite
```

### 3. RUN
In the directory that contains `main.py`, 
```Shell
(ycudb) $ python main.py
```
Once the run is complete, a SQLite database will be created in the `db` directory.

### The SQLite database contains three tables
1. new_samples
2. old_samples
3. mailed_samples

# NOTE
This script might get stuck halfway if it can't organize the entries in the input Excel files properly.

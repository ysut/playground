import pandas as pd

sheet_names = ['AD', 'Homo', 'CH', 'XL']

def df_to_excel(dfs: dataclass, output_xlsx) -> None:
    with pd.ExcelWriter(output_xlsx) as writer:
        for df, sheet_name in zip([dfs.AD, dfs.Hm, dfs.CH, dfs.XL], sheet_names):
            df.to_excel(writer, sheet_name=sheet_name, index=False)
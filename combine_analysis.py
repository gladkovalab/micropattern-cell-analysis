import pathlib
import polars as pl
from xlsxwriter import Workbook

def main():
    with Workbook("template_matching/summary.xlsx") as wb:
        for dirname, dirs, files in pathlib.Path("template_matching").walk():
            for file in files:
                if file.endswith(".xlsx") and file != "summary.xlsx":
                    sheet_name = dirname.name
                    if dirname.name == "denoised":
                        sheet_name = dirname.parent.name + "_denoised"
                    print(dirname / file)
                    df = pl.read_excel(dirname / file)
                    df.write_excel(
                        wb,
                        sheet_name[:31],
                        autofit = True
                    )



if __name__ == "__main__":
    main()

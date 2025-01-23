import duckdb

csv_files = [
    'ddg_info/ddg_info1.csv',
    'ddg_info/ddg_info2.csv',
    'ddg_info/ddg_info3.csv',
    'ddg_info/ddg_info4.csv',
    'ddg_info/ddg_info5.csv',
    'ddg_info/ddg_info6.csv',
    'ddg_info/ddg_info7.csv',
    'ddg_info/ddg_info8.csv',
    'ddg_info/ddg_info9a.csv',
    'ddg_info/ddg_info9b.csv',
    'ddg_info/ddg_info10.csv',
]
duckdb_con = duckdb.connect('ddg_info/ddg_info.db')

duckdb_con.execute(f"CREATE TABLE ddg_info AS SELECT * FROM read_csv_auto('{csv_files[0]}')")
for csv_file in csv_files[1:]:
    duckdb_con.execute(f"INSERT INTO ddg_info SELECT * FROM read_csv_auto('{csv_file}')")

duckdb_con.close()

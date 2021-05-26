# For file manipulation
import io
import zipfile

# For data structures and others
import pandas as pd


def write_dfs_to_in_memory_zip(d, **kwargs):
    arg_sep = kwargs['sep'] if 'sep' in kwargs else ';'
    arg_decimal = kwargs['decimal'] if 'decimal' in kwargs else ','
    
    mem_zip =io.BytesIO()
    
    with zipfile.ZipFile(mem_zip, mode='w', compression=zipfile.ZIP_DEFLATED) as zf:
        for key in d:
            s_buf = io.StringIO()
            d[key].to_csv(s_buf, sep=arg_sep, decimal=arg_decimal)
            s_buf.seek(0)
        
            zf.writestr(key, s_buf.getvalue())

    return mem_zip
    
    
def write_file_like_to_file(mem, otpt_file):
    with open(otpt_file, 'wb') as otpt:
        otpt.write(mem.getbuffer())
    
              
def extract_zip_to_memory(input_zip):
    input_zip=zipfile.ZipFile(input_zip)
    return {name: input_zip.open(name) for name in input_zip.namelist()}


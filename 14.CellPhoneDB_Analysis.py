
import cellphonedb
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
Result = cpdb_statistical_analysis_method.call(
        cpdb_file_path = "cellphonedb.zip",
        meta_file_path = "Meta_All_Subtype.txt",
        counts_file_path = "Exp_All_Subtype.txt",
        counts_data = 'hgnc_symbol',
        threshold = 0.1,
        output_path = ".../2.CellPhoneDB/Result")

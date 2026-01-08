import os
import shutil
import argparse
import sys
import logging
import warnings
import urllib3
import time
import traceback
from datetime import datetime

# Suppress urllib3 warnings
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
warnings.filterwarnings("ignore", category=urllib3.exceptions.InsecureRequestWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="IPython.core.display")

# Add SRC directory to Python path
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))

from second_pipeline import SecondPipeline
from enrichment_script import main as enrichment_main
from ppi_script import PPI_Analysis
from drug_gene_script import main as drug_gene_main


# ---------- Logging Helpers ----------

class AutoFlushFileHandler(logging.FileHandler):
    def emit(self, record):
        super().emit(record)
        self.flush()


class ConciseConsoleFilter(logging.Filter):
    def filter(self, record):
        msg = record.getMessage()
        return (
            record.levelname in ("ERROR", "WARNING") or
            msg.startswith(("[START]", "[STEP", "[SUCCESS]", "[INFO]", "[WELCOME]"))
        )


def setup_logging(log_file, debug=False):
    logging.getLogger('').handlers = []
    level = logging.DEBUG if debug else logging.INFO
    logger = logging.getLogger()
    logger.setLevel(level)

    file_handler = AutoFlushFileHandler(log_file, mode='w', encoding='utf-8')
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(file_handler)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(logging.Formatter('%(message)s'))
    console_handler.addFilter(ConciseConsoleFilter())
    logger.addHandler(console_handler)

    logging.getLogger('matplotlib').setLevel(logging.ERROR)
    logging.getLogger('urllib3').setLevel(logging.ERROR)

    run_id = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    logger.info("--------------------------------------------------")
    logger.info(f"[START] DeepRegulatoryNet pipeline initiated at {run_id}")
    logger.info("--------------------------------------------------")

    welcome_message = """
[WELCOME] DeepRegulatoryNet
===========================
 ____         ___ _       
|  _ \\ ___   / _ \\ \\      
| | | / __| | | | | |     
| |_| \\__ \\ | |_| | |     
|____/\\___/  \\___/|_|     
===========================
circRNA-miRNA-mRNA Pipeline
===========================
"""
    print(welcome_message)
    logger.info(welcome_message)


# ---------- Input Validation ----------

def validate_input_format(file_path, expected_prefix=None):
    with open(file_path, encoding='utf-8') as f:
        lines = [line.strip() for line in f if line.strip()]
    for line in lines:
        if expected_prefix and not line.startswith(expected_prefix):
            raise ValueError(f"Invalid ID format in {file_path}: '{line}'")
    return lines


# ---------- Main Pipeline ----------

def run_analysis(circ_file, mirna_file, deg_file, debug=False):
    logger = logging.getLogger()
    start_time = time.time()

    try:
        circ_ids = validate_input_format(circ_file, "hsa_circ_")
        mirna_ids = validate_input_format(mirna_file, "hsa-miR")
        deg_ids = validate_input_format(deg_file)
        logger.info("[INFO] Inputs: %d circRNAs, %d miRNAs, %d DEGs",
                    len(circ_ids), len(mirna_ids), len(deg_ids))

        model_file = os.path.join(os.path.dirname(__file__), 'Model_Files', 'calibrated_catboost_site_type_model.pkl')
        encoder_file = os.path.join(os.path.dirname(__file__), 'Model_Files', 'label_encoder.pkl')
        scaler_file = os.path.join(os.path.dirname(__file__), 'Model_Files', 'robust_scaler.pkl')

        for f in [circ_file, mirna_file, deg_file, model_file, encoder_file, scaler_file]:
            if not os.path.exists(f):
                raise FileNotFoundError(f"Missing file: {f}")

        for folder in ["temp", "output"]:
            if os.path.exists(folder):
                shutil.rmtree(folder)
            os.makedirs(folder)

        # -------- UPDATED ML PREDICTION ENGINE --------
        logger.info("--------------------------------------------------")
        logger.info("[STEP 1] Predicting circRNA–miRNA binding sites...")

        pipe = SecondPipeline(
            circ_file,
            mirna_file,
            deg_file,
            model_file,
            encoder_file,
            scaler_file,
            "temp",
            "output"
        )

        results = pipe.first_pipeline.process_all_circs()
        strong, allstrong = pipe.first_pipeline.find_strong_hits(results)
        pipe.first_pipeline.match_mirnas(results, strong, allstrong)

        total_circs = len(results)
        total_sites = sum(len(df) for df in results.values())
        logger.info("[INFO] Processed %d circRNAs", total_circs)
        logger.info("[INFO] Strong/Medium binding sites: %d", total_sites)

        overlap_df = pipe.first_pipeline.analyze_mrna_overlap()
        if overlap_df is not None and not overlap_df.empty:
            pipe.create_comprehensive_excel()
            logger.info("[INFO] Comprehensive interaction Excel generated")

            # Construct regulatory network (circRNA → miRNA → mRNA)
            logger.info("--------------------------------------------------")
            logger.info("[STEP 2.5] Constructing circRNA–miRNA–mRNA regulatory network...")
            try:
                G = pipe.first_pipeline.construct_network(results, strong, allstrong)
                if G is None:
                    logger.warning("[WARN] ⚠ Network construction returned no graph (missing files?)")
                else:
                    logger.info("[INFO] Regulatory network constructed: nodes=%d edges=%d", G.number_of_nodes(), G.number_of_edges())
            except Exception as e:
                logger.warning("[WARN] Network construction failed: %s", e)
                if debug:
                    logger.debug(traceback.format_exc())
        else:
            logger.warning("[WARN] ⚠ No overlapping genes found")

        # -------- DOWNSTREAM PIPELINE --------

        # Ensure overlapping genes CSV is created and get gene list
        try:
            genes = pipe.extract_overlapping_genes() or []
        except Exception as e:
            logger.warning("[WARN] Failed to extract overlapping genes: %s", e)
            if debug:
                logger.debug(traceback.format_exc())
            genes = []

        overlapping_path = os.path.join("output", "overlapping_genes.csv")

        logger.info("--------------------------------------------------")
        logger.info("[STEP 3] Performing enrichment analysis...")
        if os.path.exists(overlapping_path):
            enrichment_main(overlapping_path, output_dir=os.path.join("output", "enrichment_results"))
        else:
            logger.warning("[WARN] ⚠ Skipped enrichment: no overlapping genes file found (%s)", overlapping_path)

        logger.info("--------------------------------------------------")
        logger.info("[STEP 4] Building PPI network...")
        if os.path.exists(overlapping_path):
            PPI_Analysis(overlapping_path)
        else:
            logger.warning("[WARN] ⚠ Skipped PPI analysis: no overlapping genes file found (%s)", overlapping_path)

        logger.info("--------------------------------------------------")
        logger.info("[STEP 5] Analyzing drug–gene interactions...")
        hub_genes_path = os.path.join("output", "hub_genes.csv")
        if os.path.exists(hub_genes_path):
            drug_gene_main(hub_genes_path)
        else:
            logger.warning("[WARN] ⚠ Skipped drug–gene analysis: no hub genes file found (%s)", hub_genes_path)

        runtime = (time.time() - start_time) / 60

        logger.info("--------------------------------------------------")
        logger.info(
            f"[SUCCESS] Completed | circRNAs: {total_circs} | Sites: {total_sites} | "
            f"Overlapping genes: {len(genes)} | Runtime: {runtime:.2f} min"
        )

    except Exception as e:
        logger.error("[ERROR] Pipeline failed: %s", e)
        if debug:
            logger.error(traceback.format_exc())
        sys.exit(1)


# ---------- CLI ----------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="DeepRegulatoryNet: circRNA–miRNA–mRNA analysis pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Example usage:\n"
            "  python DeepRegulatoryNet.py --circ Test_Data/DEcircRNA.txt "
            "--mirna Test_Data/DEmiRNA.txt --deg Test_Data/DEG.txt\n\n"
            "Notes:\n"
            "  - Ensure Model_Files contains trained model and encoders.\n"
            "  - Use --debug for verbose logging and stack traces."
        )
    )
    parser.add_argument("--circ", required=True, help="Path to circRNA IDs (one per line, starting with 'hsa_circ_').")
    parser.add_argument("--mirna", required=True, help="Path to miRNA IDs (one per line, starting with 'hsa-miR').")
    parser.add_argument("--deg", required=True, help="Path to DEG gene symbols (one per line).")
    parser.add_argument("--debug", action="store_true", help="Enable debug logging and stack traces.")
    args = parser.parse_args()

    setup_logging("pipeline.log", args.debug)
    run_analysis(args.circ, args.mirna, args.deg, args.debug)

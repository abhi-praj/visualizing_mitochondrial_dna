from flask import Flask, request, send_file, jsonify
import os
from werkzeug.utils import secure_filename
from GC_Analysis import display_gcc
from SequenceLengthDistribution import sequence_length_distribution
from CodonUsageAnalysis import analyze_codon_usage
from ConservedRegions import analyze_conserved_regions
from GeneAnnotation import gene_annotation_pipeline
from PhylogeneticTreeConstruction import phylogenetic_tree_pipeline

app = Flask(__name__)
UPLOAD_FOLDER = 'uploads'
RESULTS_FOLDER = 'results'

os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(RESULTS_FOLDER, exist_ok=True)

ALLOWED_EXTENSIONS = {'fasta'}


def allowed_file(filename):
    """Check if uploaded file has an allowed extension."""
    return '.' in filename and filename.rsplit('.', 1)[
        1].lower() in ALLOWED_EXTENSIONS


@app.route('/')
def index():
    """Home route to welcome the user."""
    return "Welcome to the DNA Analysis Tool"


@app.route('/upload', methods=['POST'])
def upload_file():
    """Handle file upload and run analyses."""
    if 'file' not in request.files:
        return jsonify({'error': 'No file part'}), 400

    file = request.files['file']
    if file.filename == '':
        return jsonify({'error': 'No selected file'}), 400

    if file and allowed_file(file.filename):
        filename = secure_filename(file.filename)
        filepath = os.path.join(UPLOAD_FOLDER, filename)
        file.save(filepath)

        # Call all analysis functions here
        results = run_all_analyses(filepath)

        return jsonify(results)

    return jsonify({'error': 'File not allowed'}), 400


def run_all_analyses(filepath):
    """Run all DNA analyses and save results as image files."""
    results = {}

    # GC Content Analysis
    gc_image_path = os.path.join(RESULTS_FOLDER, 'gc_content.png')
    display_gcc(filepath, gc_image_path)
    results['gc_content'] = gc_image_path

    # Sequence Length Distribution Analysis
    seq_length_image_path = os.path.join(RESULTS_FOLDER,
                                         'sequence_length_distribution.png')
    sequence_length_distribution(filepath, seq_length_image_path)
    results['sequence_length_distribution'] = seq_length_image_path

    # Codon Usage Analysis
    codon_usage_image_path = os.path.join(RESULTS_FOLDER, 'codon_usage.png')
    analyze_codon_usage(filepath, codon_usage_image_path)
    results['codon_usage'] = codon_usage_image_path

    # Conserved Regions Analysis
    conserved_regions_image_path = os.path.join(RESULTS_FOLDER,
                                                'conserved_regions.png')
    analyze_conserved_regions(filepath, conserved_regions_image_path)
    results['conserved_regions'] = conserved_regions_image_path

    # Gene Annotation Analysis
    gene_annotation_image_folder = os.path.join(RESULTS_FOLDER,
                                                'gene_annotations')
    os.makedirs(gene_annotation_image_folder, exist_ok=True)
    gene_annotation_pipeline(filepath, gene_annotation_image_folder)
    results['gene_annotations'] = [os.path.join(gene_annotation_image_folder, f)
                                   for f in
                                   os.listdir(gene_annotation_image_folder) if
                                   f.endswith('.png')]

    # Phylogenetic Tree Construction
    phylogenetic_tree_image_path = os.path.join(RESULTS_FOLDER,
                                                'phylogenetic_tree.png')
    phylogenetic_tree_pipeline(filepath, phylogenetic_tree_image_path)
    results['phylogenetic_tree'] = phylogenetic_tree_image_path

    return results


@app.route('/download/<filename>')
def download_file(filename):
    """Allow users to download analysis result images."""
    return send_file(os.path.join(RESULTS_FOLDER, filename), as_attachment=True)


if __name__ == '__main__':
    app.run(debug=True)

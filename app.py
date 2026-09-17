#uses flask cuz yeah: picks two proteins from a already made list, will see their sequence
#identity, structural RMSD, and each one's 3D structure side by side.


from flask import Flask, render_template, jsonify, request, send_file

import alphafold_client
import sequence_compare
import structure_compare

app = Flask(__name__)

#the pre-made protein list, in the future i want to make some form of neural network 
#that pulls from alphafold but like i dunno how to do that yet whoops
PROTEINS = {
    "P69905": {
        "name": "Hemoglobin subunit alpha",
        "function": (
            "Carries oxygen from the lungs to tissues throughout the body, "
            "as part of the hemoglobin complex in red blood cells."
        ),
    },
    "P02144": {
        "name": "Myoglobin",
        "function": (
            "Stores oxygen within muscle tissue, releasing it during exertion. "
            "Evolutionarily related to hemoglobin -- both belong to the globin "
            "family and share a similar 3D fold despite serving different roles."
        ),
    },
    "P69892": {
        "name": "Hemoglobin subunit gamma-1 (fetal)",
        "function": (
            "The fetal form of hemoglobin, with a higher oxygen affinity than "
            "adult hemoglobin so the fetus can draw oxygen across the placenta."
        ),
    },
    "P01308": {
        "name": "Insulin",
        "function": (
            "A hormone that regulates blood glucose by signalling cells to "
            "take up sugar from the bloodstream."
        ),
    },
    "P61626": {
        "name": "Lysozyme C",
        "function": (
            "An enzyme found in tears, saliva and mucus that breaks down "
            "bacterial cell walls as part of the innate immune system."
        ),
    },
    "P04637": {
        "name": "Cellular tumor antigen p53",
        "function": (
            "A tumor-suppressor protein that halts cell division or triggers "
            "cell death in response to DNA damage, preventing damaged cells "
            "from proliferating."
        ),
    },
}


@app.route("/")
def index():
    return render_template("index.html", proteins=PROTEINS)


@app.route("/api/compare")
def compare():
    id_a = request.args.get("protein_a")
    id_b = request.args.get("protein_b")

    if id_a not in PROTEINS or id_b not in PROTEINS:
        return jsonify({"error": "Unknown protein ID"}), 400

    try:
        seq_a = alphafold_client.get_sequence(id_a)
        seq_b = alphafold_client.get_sequence(id_b)
        pdb_path_a = alphafold_client.get_pdb_path(id_a)
        pdb_path_b = alphafold_client.get_pdb_path(id_b)

        alignment = sequence_compare.align_sequences(seq_a, seq_b)
        identity = sequence_compare.percent_identity(alignment)

        structural_result = structure_compare.compare_structures(
            pdb_path_a, seq_a, pdb_path_b, seq_b
        )
    except (alphafold_client.AlphaFoldError, ValueError) as e:
        return jsonify({"error": str(e)}), 502

    return jsonify({
        "protein_a": {"id": id_a, **PROTEINS[id_a]},
        "protein_b": {"id": id_b, **PROTEINS[id_b]},
        "sequence_identity_percent": round(identity, 1),
        "structural_rmsd_angstrom": structural_result["rmsd_angstrom"],
        "compared_residue_pairs": structural_result["compared_residue_pairs"],
    })


@app.route("/api/structure/<uniprot_id>")
def structure_file(uniprot_id):
    #to show 3d structure on screen
    if uniprot_id not in PROTEINS:
        return "Unknown protein", 404
    try:
        path = alphafold_client.get_pdb_path(uniprot_id)
    except alphafold_client.AlphaFoldError as e:
        return str(e), 502
    return send_file(path, mimetype="chemical/x-pdb")


if __name__ == "__main__":
    app.run(debug=True, port=5000)
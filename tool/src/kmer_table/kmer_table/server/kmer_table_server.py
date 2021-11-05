from re import VERBOSE
from flask import Flask, request, jsonify, make_response
import os, sys

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import kmer_table

app = Flask(__name__, static_url_path='/')
app.config.from_pyfile('config.py')


@app.route('/')
def index():
    return app.send_static_file('index.html')

@app.route('/json')
def json():
    data = {
        "name": "ryo",
        "age": 29
    }
    return jsonify(data)

@app.route('/table_info')
def table_info():
    data = {
        "table_list": kmer_table.get_table_path_list(app.config["TABLES_DIR_PATH"]),
        "tables_dir_path" : app.config["TABLES_DIR_PATH"]
    }
    return jsonify(data)

@app.route('/gff_info')
def gff_info():
    data = {
        "gff_list": kmer_table.get_gff_path_list(app.config["GFF_DIR_PATH"]),
        "gff_dir_path" : app.config["GFF_DIR_PATH"]
    }
    return jsonify(data)

@app.route('/lookup')
def lookup():
    table_dir = request.args.get("table_path")
    target_kmer = request.args.get("kmer")
    gff_file = request.args.get("gff_path")
    feature_type = request.args.get("feature_type", default=None)
    try:
        feature_list = kmer_table.output_kmer_coord_on_feature(table_dir, kmer_table.read_gff(gff_file), len(target_kmer), target_kmer, feature_type)
    except:
        return f'"{target_kmer}" is not exist in table: {table_dir}.'
    plain_text = ""
    for f in feature_list:
        id = f["id"]
        parent = f["parent"]
        try:
            kmer_pos_list = ','.join(map(str, f['kmer_'+target_kmer]))
        except KeyError:
            continue
        plain_text += f'{f["seqid"]}\t{f["source"]}\t{f["type"]}\t{f["start"]}\t{f["end"]}\t{f["score"]}\t{f["strand"]}\t{f["phase"]}\t{f"ID={id}" if id else f"Parent={parent}"};pos_{target_kmer}={kmer_pos_list}\n'
    if plain_text == "":
        plain_text = "K-mer count 0"
    response = make_response(plain_text)
    response.mimetype = "text/plain"
    return response

@app.route('/pickup')
def pickup():
    table_dir = request.args.get("table_path")
    target_kmer = request.args.get("kmer")
    seqid_keyword = request.args.get("seqid", default="*")
    if seqid_keyword == "":
        seqid_keyword = "*"
    try:
        picked_kmer_coord_table = kmer_table.output_kmer_coord_in_table(table_dir, target_kmer,seqid_keyword)
    except:
        return f'"{target_kmer}" is not exist in table: {table_dir}.'
    plain_text = ""
    for seqid, coord_list in picked_kmer_coord_table.items():
        plain_text += f'{seqid},{",".join(map(str, coord_list))}\n'
    if plain_text == "":
        plain_text = "K-mer count 0"
    response = make_response(plain_text)
    response.mimetype = "text/plain"
    return response


if __name__ == "__main__":
    app.run(host=app.config["HOST"], port=app.config["PORT"])
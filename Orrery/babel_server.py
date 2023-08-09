from flask import Flask, request
from openbabel import openbabel

app = Flask(__name__)

@app.route('/smi_to_mol2', methods=['POST'])
def smi_to_mol2():
    data = request.get_json()
    smiles_string = data.get('smiles')

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("smi", "mol2")

    mol = openbabel.OBMol()
    obConversion.ReadString(mol, smiles_string)
    mol2_string = obConversion.WriteString(mol)

    return mol2_string

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)

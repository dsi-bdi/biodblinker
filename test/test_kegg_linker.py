import pytest
import biodblinker


@pytest.fixture(scope='module')
def linker():
    print('Setup linker')
    yield biodblinker.KEGGLinker()
    print('Teardown linker')


gene_uniprot = [
    (['mmu:69568'], [['Q6TEK5', 'Q0VGU5']]),
    (['cel:CELE_T25B9.5'], [[]]),
    (['dre:796659'], [['E7F2F5']]),
    (['spo:SPCC417.16'], [[]]),
    (['pon:100174417'], [['Q5RBB9']]),
    (['hsa:9425'], [['Q9Y232']]),
    (['ssc:397547'], [['O02835']]),
    (['xla:447703'], [['Q641H2']]),
]
@pytest.mark.parametrize("geneid,expected", gene_uniprot)
def test_convert_geneid_to_uniprot(geneid, expected, linker):
    result = linker.convert_geneid_to_uniprot(geneid)
    assert sorted(result[0]) == sorted(expected[0])


drug_drugbank = [
    (['D01300'], [['DB09006']]),
    (['D00563'], [['DB00370']]),
    (['D05209'], [['DB00957']]),
    (['D01320'], [['DB09282']]),
    (['D00345'], [['DB00808']]), 
]
@pytest.mark.parametrize("drug,expected", drug_drugbank)
def test_convert_drugid_to_drugbank(drug, expected, linker):
    result = linker.convert_drugid_to_drugbank(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_sider = [
    (['D02689'], [['CID100208902']]),
    (['D07918'], [['CID100000450']]),
    (['D09893'], [['CID111531537']]),
    (['D07866'], [['CID100003143']]),
    (['D00235'], [['CID100002249']]),
]
@pytest.mark.parametrize("drug,expected", drug_sider)
def test_convert_drugid_to_sider(drug, expected, linker):
    result = linker.convert_drugid_to_sider(drug)
    assert sorted(result[0]) == sorted(expected[0])


disease_omim = [
    (['H00246'], [['617343', '145000', '145001', '239200']]),
    (['H01900'], [['617086', '614388']]),
    (['H00976'], [['190900', '303700', '303800', '303900']]),
    (['H02386'], [[]]),
    (['H00477'], [['177170']]),
]
@pytest.mark.parametrize("disease,expected", disease_omim)
def test_convert_disease_to_omim(disease, expected, linker):
    result = linker.convert_disease_to_omim(disease)
    assert sorted(result[0]) == sorted(expected[0])


disease_mesh = [
    (['H01457'], [['D003930']]),
    (['H00174'], [['C537361', 'C565390', 'C537360']]),
    (['H01636'], [['D005356']]),
    (['H00507'], [['D019871']]),
    (['H01566'], [['D001602']]),
]
@pytest.mark.parametrize("disease,expected", disease_mesh)
def test_convert_disease_to_mesh(disease, expected, linker):
    result = linker.convert_disease_to_mesh(disease)
    assert sorted(result[0]) == sorted(expected[0])


gene_ensembl = [
    (['dre:558977'], [['ENSDARG00000069619']]),
    (['hsa:5195'], [['ENSG00000142655']]),
    (['ssc:100736663'], [['ENSSSCG00000013832']]),
    (['mmu:67830'], [['ENSMUSG00000029048']]),
    (['gga:416036'], [['ENSGALG00000006325']]),
]
@pytest.mark.parametrize("geneid,expected", gene_ensembl)
def test_convert_geneid_to_ensembl(geneid, expected, linker):
    result = linker.convert_geneid_to_ensembl(geneid)
    assert sorted(result[0]) == sorted(expected[0])


gene_names = [
    (['ath:AT5G10100'], [['TPPI']]),
    (['bta:511439'], [['SLC28A1']]),
    (['rno:689060'], [['Vom2r71']]),
    (['rno:60384'], [['p102', 'Rack2', 'Copb2']]),
    (['ath:AT3G27530'], [['GC6']]),
]
@pytest.mark.parametrize("geneid,expected", gene_names)
def test_convert_geneid_to_names(geneid, expected, linker):
    result = linker.convert_geneid_to_names(geneid)
    assert sorted(result[0]) == sorted(expected[0])


drug_chebi = [
    (['D01349'], [['32012']]),
    (['D03858'], [['4638']]),
    (['D08120'], [['63598']]),
    (['D00401'], [['50673']]),
    (['D00221'], [['28939']]),
]
@pytest.mark.parametrize("drug,expected", drug_chebi)
def test_convert_drugid_to_chebi(drug, expected, linker):
    result = linker.convert_drugid_to_chebi(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_chembl = [
    (['D06633'], [['CHEMBL2104980']]),
    (['D01128'], [['CHEMBL361812']]),
    (['D01670'], [['CHEMBL3989553']]),
    (['D01521'], [['CHEMBL3187518']]),
    (['D02829'], [['CHEMBL1200885']]),
]
@pytest.mark.parametrize("drug,expected", drug_chembl)
def test_convert_drugid_to_chembl(drug, expected, linker):
    result = linker.convert_drugid_to_chembl(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_hmdb = [
    (['D02324'], [['HMDB29940']]),
    (['D02982'], [['HMDB00517']]),
    (['D00447'], [['HMDB15150']]),
    (['D03556'], [['HMDB15561']]),
    (['D05727'], [['HMDB15338']]),
]
@pytest.mark.parametrize("drug,expected", drug_hmdb)
def test_convert_drugid_to_hmdb(drug, expected, linker):
    result = linker.convert_drugid_to_hmdb(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_hsdb = [
    (['D06114'], [['THIRAM']]),
    (['D00197'], [['RESERPINE']]),
    (['D09000'], [['TRISODIUM+PHOSPHATE']]),
    (['D00110'], [['COCAINE']]),
    (['D07958'], [['FEXOFENADINE']]),
]
@pytest.mark.parametrize("drug,expected", drug_hsdb)
def test_convert_drugid_to_hsdb(drug, expected, linker):
    result = linker.convert_drugid_to_hsdb(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_knapsack = [
    (['D00117'], [['C00029459', 'C00002654']]),
    (['D07897'], [['C00018820']]),
    (['D02798'], [['C00019560']]),
    (['D06001'], [['C00044810']]),
    (['D01035'], [['C00001836']]),
]
@pytest.mark.parametrize("drug,expected", drug_knapsack)
def test_convert_drugid_to_knapsack(drug, expected, linker):
    result = linker.convert_drugid_to_knapsack(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_ligandbox = [
    (['D09019'], [['D09019']]),
    (['D00325'], [['D00325']]),
    (['D01274'], [['D01274']]),
    (['D05139'], [['D05139']]),
    (['D01531'], [['D01531']]),
]
@pytest.mark.parametrize("drug,expected", drug_ligandbox)
def test_convert_drugid_to_ligandbox(drug, expected, linker):
    result = linker.convert_drugid_to_ligandbox(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_massbank = [
    (['D07854'], [['WA000519', 'WA000517', 'WA000516', 'WA000518', 'WA000520']]),
    (['D07064'], [['WA001214', 'WA001215', 'WA001217', 'WA001213', 'WA001218', 'WA001216']]),
    (['D07326'], [['WA001264', 'WA001267', 'WA001268', 'WA001265', 'WA001266']]),
    (['D00771'], [['WA002308', 'WA002305', 'WA002306', 'JP003097', 'WA002307']]),
    (['D08068'], [['WA001453', 'WA001449', 'WA001454', 'WA001452', 'WA001451', 'WA001450']]),
]
@pytest.mark.parametrize("drug,expected", drug_massbank)
def test_convert_drugid_to_massbank(drug, expected, linker):
    result = linker.convert_drugid_to_massbank(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_names = [
    (['D10819'], [['Aflibercept beta (genetical recombination)']]),
    (['D03602'], [['Crisnatol mesylate']]),
    (['D07762'], [['Cyhalothrin']]),
    (['D04191'], [['Fletazepam']]),
    (['D02467'], [['Droxacin sodium']]),
]
@pytest.mark.parametrize("drug,expected", drug_names)
def test_convert_drugid_to_names(drug, expected, linker):
    result = linker.convert_drugid_to_names(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_nikkaji = [
    (['D02299'], [['J17.146H']]),
    (['D01347'], [['J220.204B']]),
    (['D03027'], [['J244.505K']]),
    (['D00970'], [['J17.059C']]),
    (['D03880'], [['J2.205.416C']]),
]
@pytest.mark.parametrize("drug,expected", drug_nikkaji)
def test_convert_drugid_to_nikkaji(drug, expected, linker):
    result = linker.convert_drugid_to_nikkaji(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_pdb_ccd = [
    (['D00015'], [['GLN']]),
    (['D00204'], [['AIC']]),
    (['D03710'], [['5CH']]),
    (['D02323'], [['TOL']]),
    (['D02441'], [['EZL']]),
]
@pytest.mark.parametrize("drug,expected", drug_pdb_ccd)
def test_convert_drugid_to_pdb_ccd(drug, expected, linker):
    result = linker.convert_drugid_to_pdb_ccd(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_pubchem = [
    (['D05448'], [['47207117']]),
    (['D04324'], [['17398052']]),
    (['D02995'], [['17397151']]),
    (['D00423'], [['7847489']]),
    (['D11215'], [['384585193']]),
]
@pytest.mark.parametrize("drug,expected", drug_pubchem)
def test_convert_drugid_to_pubchem_substance(drug, expected, linker):
    result = linker.convert_drugid_to_pubchem_substance(drug)
    assert sorted(result[0]) == sorted(expected[0])


disease_names = [
    (['H02038'], [['X-linked panhypopituitarism']]),
    (['H00711'], [['Russell-Silver syndrome']]),
    (['H01138'], [['Hymenolepiasis']]),
    (['H00346'], [['Extrinsic allergic alveolitis']]),
    (['H01886'], [['Van den Ende-Gupta syndrome']]),
]
@pytest.mark.parametrize("disease,expected", disease_names)
def test_convert_disease_to_names(disease, expected, linker):
    result = linker.convert_disease_to_names(disease)
    assert sorted(result[0]) == sorted(expected[0])


disease_icd10 = [
    (['H00151'], [['E75.5']]),
    (['H00908'], [['Q87.0']]),
    (['H00035'], [['C41']]),
    (['H02294'], [['Q87.3']]),
    (['H01880'], [['Q87.8']]),
]
@pytest.mark.parametrize("disease,expected", disease_icd10)
def test_convert_disease_to_icd_10(disease, expected, linker):
    result = linker.convert_disease_to_icd_10(disease)
    assert sorted(result[0]) == sorted(expected[0])


disease_icd11 = [
    (['H00555'], [['LD2F.1Y']]),
    (['H02404'], [['1F63']]),
    (['H01502'], [['4A43.2']]),
    (['H00958'], [['9A70.Y']]),
    (['H00872'], [['LD26.4Y']]),
]
@pytest.mark.parametrize("disease,expected", disease_icd11)
def test_convert_disease_to_icd_11(disease, expected, linker):
    result = linker.convert_disease_to_icd_11(disease)
    assert sorted(result[0]) == sorted(expected[0])


disease_pubmed = [
    (['H01756'], [['21248745', '16740155', '21189955']]),
    (['H01662'], [['24326166', '17174708']]),
    (['H01344'], [['22373003', '19409520']]),
    (['H02259'], [['24619930', '24591628']]),
    (['H00277'], [['17890332', '23133618']]),
]
@pytest.mark.parametrize("disease,expected", disease_pubmed)
def test_convert_disease_to_pubmed(disease, expected, linker):
    result = linker.convert_disease_to_pubmed(disease)
    assert sorted(result[0]) == sorted(expected[0])


glycan_names = [
    (['G05525'], [['(Gal)1 (Glc)5 (GlcNAc)1 (Kdo)2']]),
    (['G04700'], [['(Gal)3 (GalNAc)1 (GlcNAc)2 (LFuc)1 (Neu5Ac)1']]),
    (['G00288'], [['LA']]),
    (['G07541'], [['(Glc)1 (Kdo)1 (Lgro-manHep)2']]),
    (['G10826'], [['(GalNAc)2 (GlcNAc)4 (Man)3 (Asn)1']]),
]
@pytest.mark.parametrize("disease,expected", glycan_names)
def test_convert_glycan_to_names(disease, expected, linker):
    result = linker.convert_glycan_to_names(disease)
    assert sorted(result[0]) == sorted(expected[0])


network_names = [
    (['N01018'], [['Mutation-caused aberrant Abeta to anterograde axonal transport']]),
    (['nt06210'], [['ERK signaling']]),
    (['N00709'], [['FAH deficiency in tyrosine degradation']]),
    (['N00892'], [['GPHN deficiency in molybdenum cofactor biosynthesis']]),
    (['N00250'], [['CDX2-overexpression to transcriptional activation']]),
]
@pytest.mark.parametrize("network,expected", network_names)
def test_convert_network_to_names(network, expected, linker):
    result = linker.convert_network_to_names(network)
    assert sorted(result[0]) == sorted(expected[0])

pathway_names = [
    (['map05340'], [['Primary immunodeficiency']]),
    (['map04745'], [['Phototransduction - fly']]),
    (['map05416'], [['Viral myocarditis']]),
    (['map00260'], [['Glycine', 'serine and threonine metabolism']]),
    (['map05221'], [['Acute myeloid leukemia']]),
]
@pytest.mark.parametrize("pathway,expected", pathway_names)
def test_convert_pathway_to_names(pathway, expected, linker):
    result = linker.convert_pathway_to_names(pathway)
    assert sorted(result[0]) == sorted(expected[0])


def test_gene_ids(linker):
    subset = {'ath:AT5G10100', 'bta:511439', 'rno:689060', 'rno:60384', 'ath:AT3G27530'}
    assert linker.gene_ids.intersection(subset) == subset


def test_drug_ids(linker):
    subset = {'D05448', 'D04324', 'D02995', 'D00423', 'D11215'}
    assert linker.drug_ids.intersection(subset) == subset


def test_disease_ids(linker):
    subset = {'H01756', 'H01662', 'H01344', 'H02259', 'H00277'}
    assert linker.disease_ids.intersection(subset) == subset


def test_network_ids(linker):
    subset = {'N01018', 'nt06210', 'N00709', 'N00892', 'N00250'}
    assert linker.network_ids.intersection(subset) == subset


def test_pathway_ids(linker):
    subset = {'map05340', 'map04745', 'map05416', 'map00260', 'map05221'}
    assert linker.pathway_ids.intersection(subset) == subset


def test_glycan(linker):
    subset = {'G05525', 'G04700', 'G00288', 'G07541', 'G10826'}
    assert linker.glycan_ids.intersection(subset) == subset

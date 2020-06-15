import pytest
import biodblinker


@pytest.fixture(scope='module')
def linker():
    print('Setup linker')
    yield biodblinker.DrugBankLinker()
    print('Teardown linker')


drug_uniprot = [
    (['DB00082'], [['P58756']]),
    (['DB04963'], [['Q05329']]),
    (['DB00030'], [['Q8HXV2']]),
    (['DB00051'], [['P01857']]),
    (['DB00070'], [['P38567']]),   
]
@pytest.mark.parametrize("drug,expected", drug_uniprot)
def test_convert_drugs_to_uniprot(drug, expected, linker):
    result = linker.convert_drugs_to_uniprot(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_kegg = [
    (['DB09216'], [['D01183']]),
    (['DB14654'], [['D01191']]),
    (['DB09014'], [['D07316']]),
    (['DB00253'], [['D02289']]),
    (['DB00184'], [['D03365']]),
]
@pytest.mark.parametrize("drug,expected", drug_kegg)
def test_convert_drugs_to_kegg_drug(drug, expected, linker):
    result = linker.convert_drugs_to_kegg_drug(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_sider = [
    (['DB06700'], [['CID100125017']]),
    (['DB01186'], [['CID100004745']]),
    (['DB00227'], [['CID100003962']]),
    (['DB00529'], [['CID100003414']]),
    (['DB06616'], [['CID105328940']]),
]
@pytest.mark.parametrize("drug,expected", drug_sider)
def test_convert_drugs_to_sider(drug, expected, linker):
    result = linker.convert_drugs_to_sider(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_bingingdb = [
    (['DB08064'], [['13345']]),
    (['DB00379'], [['50117271']]),
    (['DB09283'], [['50240032']]),
    (['DB07688'], [['11448']]),
    (['DB02379'], [['50240803']]),
]
@pytest.mark.parametrize("drug,expected", drug_bingingdb)
def test_convert_drugs_to_bindingdb(drug, expected, linker):
    result = linker.convert_drugs_to_bindingdb(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_chebi = [
    (['DB03902'], [['16995']]),
    (['DB13755'], [['31702']]),
    (['DB03795'], [['17094']]),
    (['DB01182'], [['63619']]),
    (['DB02545'], [['80003']]),
]
@pytest.mark.parametrize("drug,expected", drug_chebi)
def test_convert_drugs_to_chebi(drug, expected, linker):
    result = linker.convert_drugs_to_chebi(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_chembl = [
    (['DB15639'], [['CHEMBL3785909']]),
    (['DB08215'], [['CHEMBL1234490']]),
    (['DB12416'], [['CHEMBL1276678']]),
    (['DB11456'], [['CHEMBL1871418']]),
    (['DB13356'], [['CHEMBL2106968']]),
]
@pytest.mark.parametrize("drug,expected", drug_chembl)
def test_convert_drugs_to_chembl(drug, expected, linker):
    result = linker.convert_drugs_to_chembl(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_chemspider = [
    (['DB02377'], [['744']]),
    (['DB09395'], [['29105']]),
    (['DB03208'], [['393842']]),
    (['DB15132'], [['34980772']]),
    (['DB12762'], [['13176386']]),
]
@pytest.mark.parametrize("drug,expected", drug_chemspider)
def test_convert_drugs_to_chemspider(drug, expected, linker):
    result = linker.convert_drugs_to_chemspider(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_dpd = [
    (['DB13148'], [['4229']]),
    (['DB11342'], [['685']]),
    (['DB09561'], [['4771', '9735']]),
    (['DB09134'], [['7218', '1187']]),
    (['DB00048'], [['21378', '2150']]),
]
@pytest.mark.parametrize("drug,expected", drug_dpd)
def test_convert_drugs_to_dpd(drug, expected, linker):
    result = linker.convert_drugs_to_dpd(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_genbank = [
    (['DB00041'], [['M11144']]),
    (['DB00002'], [['J00228']]),
    (['DB00032'], [['X00264']]),
    (['DB00034'], [['J00207']]),
    (['DB00052'], [['AF374232']]),
]
@pytest.mark.parametrize("drug,expected", drug_genbank)
def test_convert_drugs_to_genbank(drug, expected, linker):
    result = linker.convert_drugs_to_genbank(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_pharma = [
    (['DB00734'], [['96']]),
    (['DB00747'], [['330']]),
    (['DB00841'], [['535']]),
    (['DB00508'], [['4330']]),
    (['DB04466'], [['2763']]),
]
@pytest.mark.parametrize("drug,expected", drug_pharma)
def test_convert_drugs_to_guide_pharmacology(drug, expected, linker):
    result = linker.convert_drugs_to_guide_pharmacology(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_iuphar = [
    (['DB01182'], [['2561']]),
    (['DB01085'], [['305']]),
    (['DB00657'], [['3990']]),
    (['DB00692'], [['502']]),
    (['DB01054'], [['2334']]),
]
@pytest.mark.parametrize("drug,expected", drug_iuphar)
def test_convert_drugs_to_iuphar(drug, expected, linker):
    result = linker.convert_drugs_to_iuphar(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_pdb = [
    (['DB03295'], [['TS5']]),
    (['DB03248'], [['PEZ']]),
    (['DB02258'], [['254']]),
    (['DB08632'], [['TMM']]),
    (['DB03668'], [['BMQ']]),
]
@pytest.mark.parametrize("drug,expected", drug_pdb)
def test_convert_drugs_to_pdb(drug, expected, linker):
    result = linker.convert_drugs_to_pdb(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_pharmgkb = [
    (['DB00740'], [['PA451251']]),
    (['DB01024'], [['PA164748728']]),
    (['DB00054'], [['PA448006']]),
    (['DB00202'], [['PA451522']]),
    (['DB01584'], [['PA164754910']]),
]
@pytest.mark.parametrize("drug,expected", drug_pharmgkb)
def test_convert_drugs_to_pharmgkb(drug, expected, linker):
    result = linker.convert_drugs_to_pharmgkb(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_pubchem_cmpd = [
    (['DB04850'], [['213049']]),
    (['DB03882'], [['242332']]),
    (['DB02983'], [['65258']]),
    (['DB00441'], [['60750']]),
    (['DB03608'], [['2354']]),
]
@pytest.mark.parametrize("drug,expected", drug_pubchem_cmpd)
def test_convert_drugs_to_pubchem_compound(drug, expected, linker):
    result = linker.convert_drugs_to_pubchem_compound(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_pubchem_sub = [
    (['DB03129'], [['46506340']]),
    (['DB13136'], [['347829254']]),
    (['DB08226'], [['99444697']]),
    (['DB03734'], [['46508938']]),
    (['DB04410'], [['46506901']]),
]
@pytest.mark.parametrize("drug,expected", drug_pubchem_sub)
def test_convert_drugs_to_pubchem_substance(drug, expected, linker):
    result = linker.convert_drugs_to_pubchem_substance(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_ttd = [
    (['DB00274'], [['DAP001179']]),
    (['DB01047'], [['DAP000421']]),
    (['DB00047'], [['DAP001088']]),
    (['DB01067'], [['DAP000920']]),
    (['DB04896'], [['DAP001155']]),
]
@pytest.mark.parametrize("drug,expected", drug_ttd)
def test_convert_drugs_to_ttd(drug, expected, linker):
    result = linker.convert_drugs_to_ttd(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_wikipedia = [
    (['DB14027'], [['Taspoglutide']]),
    (['DB00568'], [['Cinnarizine']]),
    (['DB00359'], [['Sulfadiazine']]),
    (['DB02045'], [['Pentylamine']]),
    (['DB02057'], [['AMPA']]),
]
@pytest.mark.parametrize("drug,expected", drug_wikipedia)
def test_convert_drugs_to_wikipedia(drug, expected, linker):
    result = linker.convert_drugs_to_wikipedia(drug)
    assert sorted(result[0]) == sorted(expected[0])

drug_sec_id = [
    (['DB01887'], [['DB01887', 'EXPT00659']]),
    (['DB03018'], [['DB03018', 'EXPT00155']]),
    (['DB00641'], [['DB00641', 'APRD00104']]),
    (['DB03167'], [['EXPT02511', 'DB03167']]),
    (['DB02579'], [['DB02579', 'EXPT00484']]),
]
@pytest.mark.parametrize("drug,expected", drug_sec_id)
def test_convert_drugs_to_secondaryids(drug, expected, linker):
    result = linker.convert_drugs_to_secondaryids(drug)
    assert sorted(result[0]) == sorted(expected[0])


sec_id_drug = [
    (['EXPT02173'], [['DB02467']]),
    (['APRD00061'], [['DB01176']]),
    (['APRD00789'], [['DB01163']]),
    (['EXPT01276'], [['DB03031']]),
    (['EXPT00104'], [['DB04313']]),
]
@pytest.mark.parametrize("sec_id,expected", sec_id_drug)
def test_convert_secondaryids_to_drugs(sec_id, expected, linker):
    result = linker.convert_secondaryids_to_drugs(sec_id)
    assert sorted(result[0]) == sorted(expected[0])

name_drugs = [
    (['n_gamma__nitro_l_arginine'], [['DB04223']]),
    (['ertugliflozin'], [['DB11827']]),
    (['2_amino_3__4_hydroxy_6_oxo_3__2_phenyl_cyclopropylimino__cyclohexa_1_4_dienyl__propionic_acid'], [['DB01657']]),
    (['good_health_formula'], [['DB00126']]),
    (['alfameprodina'], [['DB01499']]),
]
@pytest.mark.parametrize("name,expected", name_drugs)
def test_convert_names_to_drugs(name, expected, linker):
    result = linker.convert_names_to_drugs(name)
    assert sorted(result[0]) == sorted(expected[0])

drug_names = [
    (['DB02133'], [['Chlorophyll A']]),
    (['DB08840'], [['Nicotinic acid methylamide', 'N-Methyl-3-pyridinecarboxamide', '3-(N-Methylcarbamoyl)pyridine', '3-(Methylcarbamoyl)pyridine', 'N-methylnicotinamide', 'N-Methyl nicotineamide', 'Nicotinyl methylamide', 'N-methylpyridine-3-carboxamide']]),
    (['DB05496'], [['Otelixizumab']]),
    (['DB14750'], [['Cidoxepina', 'Cidoxepinum', '(Z)-Doxepin', 'Cidoxepin']]),
    (['DB03487'], [['(S)-Aspartimide', '3-aminosuccinimide']]),
]
@pytest.mark.parametrize("drug,expected", drug_names)
def test_convert_drugs_to_names(drug, expected, linker):
    result = linker.convert_drugs_to_names(drug)
    assert sorted(result[0]) == sorted(expected[0])


def test_drug_ids(linker):
    subset = {'DB09229', 'DB02545', 'DB07688', 'DB09196', 'DB00844'}
    assert linker.drug_ids.intersection(subset) == subset
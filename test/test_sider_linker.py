import pytest
import biodblinker


@pytest.fixture(scope='module')
def linker():
    print('Setup linker')
    yield biodblinker.SiderLinker()
    print('Teardown linker')


drug_kegg = [
    (['CID100005359'], [['D00452']]),
    (['CID100000750'], [['D00011']]),
    (['CID100060843'], [['D07472', 'D03828', 'D06503']]),
    (['CID100005379'], [['D08011', 'D00589']]),
    (['CID100003372'], [['D07977', 'D02163', 'D00791']]),
]
@pytest.mark.parametrize("drug,expected", drug_kegg)
def test_convert_drugs_to_kegg_drug(drug, expected, linker):
    result = linker.convert_drugs_to_kegg_drug(drug)
    assert sorted(result[0]) == sorted(expected[0])


drug_drugbank = [
    (['CID100000232'], [['DB04027', 'DB00125']]),
    (['CID100004197'], [['DB00235']]),
    (['CID100003367'], [['DB04441']]),
    (['CID100003382'], [['DB01047']]),
    (['CID100005726'], [['DB00495']]),
]
@pytest.mark.parametrize("drug,expected", drug_drugbank)
def test_convert_drugs_to_drugbank(drug, expected, linker):
    result = linker.convert_drugs_to_drugbank(drug)
    assert sorted(result[0]) == sorted(expected[0])

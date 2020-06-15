import pytest
import biodblinker


@pytest.fixture(scope='module')
def linker():
    print('Setup linker')
    yield biodblinker.UniprotLinker()
    print('Teardown linker')


def test_uniprot_ids(linker):
    subset = {'P10636', 'Q01525', 'Q84J55', 'P41932', 'P32234'}
    acc_ids = linker.uniprot_ids
    assert acc_ids.intersection(subset) == subset


id_names = [
    (['P10636'], [['TAU_HUMAN']]),
    (['Q01525'], [['14332_ARATH']]),
    (['Q84J55'], [['14331_ORYSJ']]),
    (['P41932'], [['14331_CAEEL']]),
    (['P32234'], [['128UP_DROME']]),
]
@pytest.mark.parametrize("acc,expected", id_names)
def test_uniprot_names(acc, expected, linker):
    result = linker.convert_uniprot_to_names(acc)
    assert sorted(result[0]) == sorted(expected[0])


id_kegg_genes = [
    (['D3ZKP1'], [['rno:499506']]),
    (['A0A1L8HIW1'], [['xla:108708675']]),
    (['Q5PQ60'], [['xla:495972']]),
    (['Q0P5D8'], [['bta:533791']]),
    (['Q3TUY3'], [['mmu:638345']]),
]
@pytest.mark.parametrize("acc,expected", id_kegg_genes)
def test_convert_uniprot_to_kegg(acc, expected, linker):
    result = linker.convert_uniprot_to_kegg(acc)
    assert sorted(result[0]) == sorted(expected[0])


id_drugbank = [
    (['P58756'], [['DB00082', 'DB00052']]),
    (['Q05329'], [['DB04963']]),
    (['Q8HXV2'], [['DB00030', 'DB00071']]),
    (['P01857'], [['DB00051', 'DB00073', 'DB00043',
     'DB00065', 'DB00074', 'DB00072', 'DB00081', 'DB00111']]),
    (['P38567'], [['DB00070']]),
]
@pytest.mark.parametrize("acc,expected", id_drugbank)
def test_convert_uniprot_to_drugbank(acc, expected, linker):
    result = linker.convert_uniprot_to_drugbank(acc)
    assert sorted(result[0]) == sorted(expected[0])


id_hpa = [
    (['K7EJ46'], [['HPA077331']]),
    (['Q9BRK5'], [['HPA011249', 'HPA011895', 'CAB015227']]),
    (['Q9UGY1'], [['HPA003547']]),
    (['Q8N0Z3'], [['HPA064843']]),
    (['Q494X1'], [['HPA046691', 'HPA052337', 'HPA042585']]),
]
@pytest.mark.parametrize("acc,expected", id_hpa)
def test_convert_uniprot_to_hpa(acc, expected, linker):
    result = linker.convert_uniprot_to_hpa(acc)
    assert sorted(result[0]) == sorted(expected[0])


id_biogrid = [
    (['Q15691'], [['116581']]),
    (['P37691'], [['4263302']]),
    (['Q5ZKF5'], [['677840']]),
    (['P53036'], [['33010']]),
    (['Q6NVH7'], [['125954']]),
]
@pytest.mark.parametrize("acc,expected", id_biogrid)
def test_convert_uniprot_to_biogrid(acc, expected, linker):
    result = linker.convert_uniprot_to_biogrid(acc)
    assert sorted(result[0]) == sorted(expected[0])


id_chembl = [
    (['P15104'], [['CHEMBL4612']]),
    (['P50637'], [['CHEMBL2149']]),
    (['P80365'], [['CHEMBL3746']]),
    (['Q8VCH6'], [['CHEMBL3774292']]),
    (['Q5KTC7'], [['CHEMBL2034805']]),
]
@pytest.mark.parametrize("acc,expected", id_chembl)
def test_convert_uniprot_to_chembl(acc, expected, linker):
    result = linker.convert_uniprot_to_chembl(acc)
    assert sorted(result[0]) == sorted(expected[0])


id_dip  = [
    (['Q12358'], [['DIP-1933N']]),
    (['Q84M91'], [['DIP-46424N']]),
    (['P28743'], [['DIP-3012N']]),
    (['Q9VCI0'], [['DIP-17679N']]),
    (['P33315'], [['DIP-757N']]),
]
@pytest.mark.parametrize("acc,expected", id_dip)
def test_convert_uniprot_to_dip(acc, expected, linker):
    result = linker.convert_uniprot_to_dip(acc)
    assert sorted(result[0]) == sorted(expected[0])


id_embl_cds = [
    (['Q8TGR1'], [['AAL79240.1']]),
    (['P0A834'], [['AAG57730.1', 'BAB36905.1']]),
    (['Q9Z0F3'], [['AAD08703.1', 'AAC83150.1', 'BAE22856.1', 'BAE26403.1', 'BAE36742.1', 'EDL26315.1', 'AAH52690.1']]),
    (['Q16790'], [['CAA47315.1', 'CAB82444.1', 'EAW58359.1', 'AAH14950.1']]),
    (['P49185'], [['AAA42111.1']]),
]
@pytest.mark.parametrize("acc,expected", id_embl_cds)
def test_convert_uniprot_to_embl_cds(acc, expected, linker):
    result = linker.convert_uniprot_to_embl_cds(acc)
    assert sorted(result[0]) == sorted(expected[0])


id_embl = [
    (['Q9DG67'], [['AF178529']]),
    (['P97857'], [['AB001735', 'D67076', 'AC126936', 'BC040382', 'BC050834']]),
    (['Q12358'], [['Z47973', 'Z73162', 'AY692737', 'BK006945']]),
    (['P9WJV5'], [['AL123456']]),
    (['Q0DBW8'], [['AP008212', 'AP014962', 'AK119453']]),
]
@pytest.mark.parametrize("acc,expected", id_embl)
def test_convert_uniprot_to_embl(acc, expected, linker):
    result = linker.convert_uniprot_to_embl(acc)
    assert sorted(result[0]) == sorted(expected[0])


id_ensembl = [
    (['P23772'], [['ENSMUSG00000015619']]),
    (['Q2KHV2'], [['ENSBTAG00000000758']]),
    (['Q2TGJ4'], [['ENSRNOG00000002751']]),
    (['P47992'], [['ENSG00000143184']]),
    (['P32443'], [['ENSMUSG00000036144']]),
]
@pytest.mark.parametrize("acc,expected", id_ensembl)
def test_convert_uniprot_to_ensembl(acc, expected, linker):
    result = linker.convert_uniprot_to_ensembl(acc)
    assert sorted(result[0]) == sorted(expected[0])


id_gene_syn = [
    (['Q3TT99'], [['Gde6']]),
    (['F1M391'], [['Tmem173']]),
    (['Q91YX0'], [['Icb1']]),
    (['P54115'], [['ALDH1']]),
    (['Q6NZA9'], [['Taf9l']]),
]
@pytest.mark.parametrize("acc,expected", id_gene_syn)
def test_convert_uniprot_to_gene_synonym(acc, expected, linker):
    result = linker.convert_uniprot_to_gene_synonym(acc)
    assert sorted(result[0]) == sorted(expected[0])


id_gi = [
    (['Q80W21'], [['74196865', '74193997', '148669982', '113679874', '81912821', '110591195', '30354091', '110591196']]),
    (['P53680'], [['119577852', '13623469', '189065218', '1296607', '666637993', '70906430', '3413477', '666637995', '3413475', '666637997', '51338780', '11038643']]),
    (['A0A023PZJ9'], [['614461535', '728050414']]),
    (['Q54N41'], [['60466546', '74853808', '66808759']]),
    (['Q27571'], [['1000082', '72151194', '190358925', '33589436', '24583543', '665407497', '599125272', '6707649', '78706870', '7297764', '26006989']]),
]
@pytest.mark.parametrize("acc,expected", id_gi)
def test_convert_uniprot_to_gi(acc, expected, linker):
    result = linker.convert_uniprot_to_gi(acc)
    assert sorted(result[0]) == sorted(expected[0])


id_hgnc = [
    (['Q14582'], [['HGNC:13906']]),
    (['Q9UQF0'], [['HGNC:13525']]),
    (['Q8NET5'], [['HGNC:29872']]),
    (['P55895'], [['HGNC:9832']]),
    (['Q8WTT2'], [['HGNC:24034']]),
]
@pytest.mark.parametrize("acc,expected", id_hgnc)
def test_convert_uniprot_to_hgnc(acc, expected, linker):
    result = linker.convert_uniprot_to_hgnc(acc)
    assert sorted(result[0]) == sorted(expected[0])


id_mim = [
    (['P56747'], [['615798']]),
    (['Q9NYW5'], [['604869']]),
    (['Q8IZL8'], [['609455']]),
    (['Q8IY26'], [['611666']]),
    (['Q8IXV7'], [['236000', '613169']]),
]
@pytest.mark.parametrize("acc,expected", id_mim)
def test_convert_uniprot_to_mim(acc, expected, linker):
    result = linker.convert_uniprot_to_mim(acc)
    assert sorted(result[0]) == sorted(expected[0])


id_mint = [
    (['Q9R059'], [['Q9R059']]),
    (['Q8N8H1'], [['Q8N8H1']]),
    (['P28274'], [['P28274']]),
    (['P00669'], [['P00669']]),
    (['Q8VZS8'], [['Q8VZS8']]),
]
@pytest.mark.parametrize("acc,expected", id_mint)
def test_convert_uniprot_to_mint(acc, expected, linker):
    result = linker.convert_uniprot_to_mint(acc)
    assert sorted(result[0]) == sorted(expected[0])


id_pharmgkb = [
    (['Q9BYW3'], [['PA27246']]),
    (['Q96PE2'], [['PA134884537']]),
    (['P52803'], [['PA27660']]),
    (['Q9NVR2'], [['PA142672354']]),
    (['Q7Z4H4'], [['PA134898869']]),
]
@pytest.mark.parametrize("acc,expected", id_pharmgkb)
def test_convert_uniprot_to_pharmgkb(acc, expected, linker):
    result = linker.convert_uniprot_to_pharmgkb(acc)
    assert sorted(result[0]) == sorted(expected[0])


id_protemoicsdb = [
    (['O95400'], [['50853']]),
    (['Q8WUY1'], [['74721']]),
    (['Q9Y4Z0'], [['86267']]),
    (['P34096'], [['54939']]),
    (['P62750'], [['57421']]),
]
@pytest.mark.parametrize("acc,expected", id_protemoicsdb)
def test_convert_uniprot_to_proteomicsdb(acc, expected, linker):
    result = linker.convert_uniprot_to_proteomicsdb(acc)
    assert sorted(result[0]) == sorted(expected[0])


id_refseq = [
    (['O08394'], [['NP_388606.1', 'WP_003242884.1']]),
    (['Q99J99'], [['NP_001155964.1', 'NP_001155965.1', 'NP_619611.3', 'XP_006521057.1', 'XP_006521058.1', 'XP_006521059.1']]),
    (['G5EGM3'], [['NP_493426.2']]),
    (['Q8NGK1'], [['NP_001005237.1']]),
    (['O43298'], [['NP_001129248.1', 'NP_054726.1', 'XP_005251892.1', 'XP_011516711.1', 'XP_011516713.1']]),
]
@pytest.mark.parametrize("acc,expected", id_refseq)
def test_convert_uniprot_to_refseq(acc, expected, linker):
    result = linker.convert_uniprot_to_refseq(acc)
    assert sorted(result[0]) == sorted(expected[0])


id_string = [
    (['P0AFW6'], [['155864.EDL933_0683']]),
    (['Q01066'], [['10116.ENSRNOP00000052147']]),
    (['P9WHW1'], [['83332.Rv2582']]),
    (['Q61884'], [['10090.ENSMUSP00000034746']]),
    (['P16043'], [['10090.ENSMUSP00000029172']]),
]
@pytest.mark.parametrize("acc,expected", id_string)
def test_convert_uniprot_to_string(acc, expected, linker):
    result = linker.convert_uniprot_to_string(acc)
    assert sorted(result[0]) == sorted(expected[0])


id_uniparc = [
    (['O94458'], [['UPI000228F407']]),
    (['Q07876'], [['UPI0000060376']]),
    (['Q9VW12'], [['UPI000007E777']]),
    (['Q9ZVI6'], [['UPI00000A3665']]),
    (['Q499U1'], [['UPI000017ED9E']]),
]
@pytest.mark.parametrize("acc,expected", id_uniparc)
def test_convert_uniprot_to_uniparc(acc, expected, linker):
    result = linker.convert_uniprot_to_uniparc(acc)
    assert sorted(result[0]) == sorted(expected[0])


id_uniprotkb = [
    (['P0A7B2'], [['PPK1_ECO57']]),
    (['Q10M50'], [['CHLH_ORYSJ']]),
    (['Q5RJY4'], [['DRS7B_RAT']]),
    (['Q8CBY0'], [['GATC_MOUSE']]),
    (['Q9SEZ3'], [['CDF5_ARATH']]),
]
@pytest.mark.parametrize("acc,expected", id_uniprotkb)
def test_convert_uniprot_to_uniprotkb_id(acc, expected, linker):
    result = linker.convert_uniprot_to_uniprotkb_id(acc)
    assert sorted(result[0]) == sorted(expected[0])


id_uniref = [
    (['P0A925'], [['UniRef100_P0A925']]),
    (['P55868'], [['UniRef100_P55868']]),
    (['P96668'], [['UniRef100_P96668']]),
    (['Q9UU82'], [['UniRef100_Q9UU82']]),
    (['Q9FK28'], [['UniRef100_Q9FK28']]),
]
@pytest.mark.parametrize("acc,expected", id_uniref)
def test_convert_uniprot_to_uniref100(acc, expected, linker):
    result = linker.convert_uniprot_to_uniref100(acc)
    assert sorted(result[0]) == sorted(expected[0])

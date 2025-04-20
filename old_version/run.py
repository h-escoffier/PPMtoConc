from old_version.mw import load_from_paxdb
from ppm_to_conc.conversion import convert_data


def test(): 
    run('input/9606-iBAQ_MCF7_Geiger_2012_uniprot.txt', 'output/MCF7_test.txt', is_paxdb=True)


def run(path_input, path_output, total_prot_content=None, is_paxdb=True):
    """
    
    """ 
    if is_paxdb: 
        # Load data from PaxDBs
        load_from_paxdb(path_input, path_output)
    # Convert data and save to CSV
    if total_prot_content is None:
        convert_data(path_output, path_output)
    else:
        convert_data(path_output, path_output, total_prot_content)


if __name__ == "__main__":
    print('start')
    test()
    print('end')
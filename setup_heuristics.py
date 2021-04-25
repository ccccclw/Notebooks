#!/usr/bin/env python
# encoding: utf-8

import numpy as np
from meld.remd import ladder, adaptor, master_runner
from meld import comm, vault
from meld import system
from meld import parse
import meld.system.montecarlo as mc
from meld.system.restraints import LinearRamp,ConstantRamp
from collections import namedtuple
import glob as glob 


N_REPLICAS = 30
N_STEPS = 10000
BLOCK_SIZE = 100


hydrophobes = 'AILMFPWV'
hydrophobes_res = ['ALA','ILE','LEU','MET','PHE','PRO','TRP','VAL']

def create_Evolution(s,accuracy=0.7,scaler=None,fname='contacts.dat'):
    scaler = scaler if scaler else s.restraints.create_scaler('nonlinear', alpha_min=0.4, alpha_max=1.0, factor=4.0)
    Evolution_restraints = []
    #Heavy atoms for each residue
    atoms = {"ALA":['N','C','O','CA','CB'],
             "VAL":['N','C','O','CA','CB','CG1','CG2'],
             "LEU":['N','C','O','CA','CB','CG','CD1','CD2'],
             "ILE":['N','C','O','CA','CB','CG1','CG2','CD1'],
             "PHE":['N','C','O','CA','CB','CG','CD1','CE1','CZ','CE2','CD2'],
             "TRP":['N','C','O','CA','CB','CG','CD1','NE1','CE2','CZ2','CH2','CZ3','CE3','CD2'],
             "MET":['N','C','O','CA','CB','CG','SD','CE'],
             "PRO":['N','C','O','CD','CG','CB','CA'],
             "ASP":['N','C','O','CA','CB','CG','OD1','OD2'],
             "GLU":['N','C','O','CA','CB','CG','CD','OE1','OE2'],
             "LYS":['N','C','O','CA','CB','CG','CD','CE','NZ'],
             "ARG":['N','C','O','CA','CB','CG','CD','NE','CZ','NH1','NH2'],
             "HIS":['N','C','O','CA','CB','CG','ND1','CE1','NE2','CD2'],
             "HID":['N','C','O','CA','CB','CG','ND1','CE1','NE2','CD2'],
             "HIE":['N','C','O','CA','CB','CG','ND1','CE1','NE2','CD2'],
             "HIP":['N','C','O','CA','CB','CG','ND1','CE1','NE2','CD2'],
             "GLY":['N','C','O','CA'],
             "SER":['N','C','O','CA','CB','OG'],
             "THR":['N','C','O','CA','CB','CG2','OG1'],
             "CYS":['N','C','O','CA','CB','SG'],
             "CYX":['N','C','O','CA','CB','SG'],
             "TYR":['N','C','O','CA','CB','CG','CD1','CE1','CZ','OH','CE2','CD2'],
             "ASN":['N','C','O','CA','CB','CG','OD1','ND2'],
             "GLN":['N','C','O','CA','CB','CG','CD','OE1','NE2']}

    #s.residue_numbers and s.residue_names have as many instances as atoms
    #we first use zip to create tuples and then set to create unique (disordered) lists of tuples
    #the sorted organizes them. Finally the dict creates an instance (one based) of residue to residue name
    sequence = [(i,j) for i,j in zip(s.residue_numbers,s.residue_names)]
    sequence = sorted(set(sequence))
    sequence = dict(sequence)

    #Read file, 1 based
    lines = open(fname).read().splitlines()
    for i,line in enumerate(lines):
        cols = line.split()
        index_i = int(cols[0])
        index_j = int(cols[1])
        res_i = sequence[index_i]
        res_j = sequence[index_j]
        atoms_i = atoms[res_i]
        atoms_j = atoms[res_j]

        local_contact = []
        for a_i in atoms_i:
            for a_j in atoms_j:
                #EvFold paper defines contact as 5angstrong any atom contact. We are only using heavy atoms
                #Hence, should allow to get closer.
                local_contact.append(s.restraints.create_restraint('distance', scaler, LinearRamp(0,100,0,1),r1=0.0, r2=0.0, r3=0.60, r4=0.75, k=250.0,
            atom_1_res_index=index_i, atom_1_name=a_i, atom_2_res_index=index_j, atom_2_name=a_j))

        Evolution_restraints.append(s.restraints.create_restraint_group(local_contact,1))

        #print 'Evolution:',index_i,index_j,res_i,res_j

    all_rest = len(Evolution_restraints)
    active = int( accuracy * all_rest )
    print(active,all_rest)
    s.restraints.add_selectively_active_collection(Evolution_restraints, active)

def choose_hydrophobic_factor(ContactOrder):
    factor = 1.8834*np.log(ContactOrder) - 2.4853
    if factor < 1.:
       factor = 1.
    return factor
    
def choose_strandpair_factor(ContactOrder):
    factor = 1.6856*np.log(ContactOrder) - 3.0749
    if factor < 1.:
       factor = 1.
    return factor
    
def choose_hbond_factor(ContactOrder):
    factor = 1.6691*np.log(ContactOrder) - 1.9437
    if factor < 1.:
       factor = 1.
    return factor

class make_CO_scaler(object):
    def __init__(self,system=None,scaler_type='plateaunonlinear',alpha_min=0,alpha_one=0.25,alpha_two=0.5,alpha_max=1.0,strength_at_alpha_min=1.0, strength_at_alpha_max=0.0):
        self._system = system
        self._scaler_type = scaler_type
        self._alpha_min = alpha_min
        self._alpha_one = alpha_one
        self._alpha_two = alpha_two
        self._alpha_max = alpha_max
        self._strength_at_alpha_min = strength_at_alpha_min
        self._strength_at_alpha_max = strength_at_alpha_max


    def __call__(self,CO,type):
        if (type is 'hydrophobic'):
           factor = choose_hydrophobic_factor(CO)
        if (type is 'strand'):
           factor = choose_strandpair_factor(CO)
        if (type is 'Hbond'):
           factor = choose_hbond_factor(CO)
        if self._alpha_one:
            #print type,self._alpha_min,self._alpha_max
            self._system.restraints.create_scaler(self._scaler_type, alpha_min=self._alpha_min, alpha_one=self._alpha_one, alpha_two=self._alpha_two, alpha_max=self._alpha_max, factor=factor, strength_at_alpha_min=self._strength_at_alpha_min, strength_at_alpha_max=self._strength_at_alpha_max)
        else:
            self._system.restraints.create_scaler(self._scaler_type, alpha_min=self._alpha_min, alpha_max=self._alpha_max, factor=factor, strength_at_alpha_min=self._strength_at_alpha_min, strength_at_alpha_max=self._strength_at_alpha_max)


def make_ss_groups(subset=None):
           active = 0
           extended = 0
           sse = []
           ss = open('ss.dat','r').readlines()[0]
           for i,l in enumerate(ss.rstrip()):
               #print i,l
               if l not in "HE.":
                   continue
               if l not in 'E' and extended:
                   end = i
                   sse.append((start+1,end))
                   extended = 0
               if l in 'E':
                   if i+1 in subset:
                       active = active + 1
                   if extended:
                       continue
                   else:
                       start = i
                       extended = 1
           #print active
           return sse,active


def _get_secondary_sequence(filename=None, contents=None, file=None):
    contents = parse._handle_arguments(filename, contents, file)
    lines = contents.splitlines()
    lines = [line.strip() for line in lines if not line.startswith('#')]
    sequence = ''.join(lines)
    for ss in sequence:
        if not ss in 'HE.':
            raise RuntimeError('Unknown secondary structure type "{}"'.format(aa))
    return sequence



def generate_strand_pairs(s,sse,active,scaler,active_per_cent=0.45,subset=np.array([]),CO=True):
    n_res = s.residue_numbers[-1]
    subset = subset if subset.size else np.array(list(range(n_res)))+1 
    strand_pair = []
    for i in range(len(sse)):
        start_i,end_i = sse[i]
        for j in range(i+1,len(sse)):
            start_j,end_j = sse[j]
            
            for res_i in range(start_i,end_i+1):
                for res_j in range(start_j,end_j+1):
                    if res_i in subset or res_j in subset:
                        #print res_i,res_j
                        g = []
                        make_pairNO(g,s,res_i,res_j,scaler,CO)
                        strand_pair.append(s.restraints.create_restraint_group(g,1))
                        g = []
                        make_pairON(g,s,res_i,res_j,scaler,CO)
                        strand_pair.append(s.restraints.create_restraint_group(g,1))
    all_rest = len(strand_pair)
    active = int(active * active_per_cent)
    print("strand_pairs:", all_rest,active)
    s.restraints.add_selectively_active_collection(strand_pair, active)

def make_pairNO(g,s,i,j,scaler,CO=True):
    contact_order = abs(i-j)
    if CO:
        tmp_scaler = scaler(contact_order,'strand')
    else:
        tmp_scaler = scaler
    g.append(s.restraints.create_restraint('distance', tmp_scaler, LinearRamp(0,100,0,1),r1=0.0, r2=0.0, r3=0.35, r4=0.50, k=250.0,
            atom_1_res_index=i, atom_1_name='N', atom_2_res_index=j, atom_2_name='O'))

def make_pairON(g,s,i,j,scaler,CO=True):
    contact_order = abs(i-j)
    if CO:
        tmp_scaler = scaler(contact_order,'strand')
    else:
        tmp_scaler = scaler
    g.append(s.restraints.create_restraint('distance', tmp_scaler,LinearRamp(0,100,0,1), r1=0.0, r2=0.0, r3=0.35, r4=0.50, k=250.0,
            atom_1_res_index=i, atom_1_name='O', atom_2_res_index=j, atom_2_name='N'))

def create_hydrophobes(s,ContactsPerHydroph=1.3,scaler=None,group_1=np.array([]),group_2=np.array([]),CO=True):
    atoms = {"ALA":['CA','CB'],
             "VAL":['CA','CB','CG1','CG2'],
             "LEU":['CA','CB','CG','CD1','CD2'],
             "ILE":['CA','CB','CG1','CG2','CD1'],
             "PHE":['CA','CB','CG','CD1','CE1','CZ','CE2','CD2'],
             "TRP":['CA','CB','CG','CD1','NE1','CE2','CZ2','CH2','CZ3','CE3','CD2'],
             "MET":['CA','CB','CG','SD','CE'],
             "PRO":['CD','CG','CB','CA']}
    #Groups should be 1 centered
    n_res = s.residue_numbers[-1]
    group_1 = group_1 if group_1.size else np.array(list(range(n_res)))+1
    group_2 = group_2 if group_2.size else np.array(list(range(n_res)))+1
    scaler = scaler if scaler else s.restraints.create_scaler('nonlinear', alpha_min=0.4, alpha_max=1.0, factor=4.0)

    #Get a list of names and residue numbers, if just use names might skip some residues that are two
    #times in a row
    #make list 1 centered
    sequence = [(i,j) for i,j in zip(s.residue_numbers,s.residue_names)]
    sequence = sorted(set(sequence))
    sequence = dict(sequence)

    #print sequence
    #print hydrophobes_res
    #Get list of groups with only residues that are hydrophobs
    group_1 = [ res for res in group_1 if (sequence[res] in hydrophobes_res) ]
    group_2 = [ res for res in group_2 if (sequence[res] in hydrophobes_res) ]


    pairs = []
    hydroph_restraints = []
    for i in group_1:
        for j in group_2:

            # don't put the same pair in more than once
            if ( (i,j) in pairs ) or ( (j,i) in pairs ):
                continue

            if ( i ==j ):
                continue

            if (abs(i-j)< 7):
                continue
            pairs.append( (i,j) )
           
            atoms_i = atoms[sequence[i]]
            atoms_j = atoms[sequence[j]]
   
            local_contact = []
            for a_i in atoms_i:
                for a_j in atoms_j:
                    if CO:
                        #print i,j,a_i,a_j
                        tmp_scaler = scaler(abs(i-j), 'hydrophobic')
                    else:
                        tmp_scaler = scaler
                    local_contact.append(s.restraints.create_restraint('distance', tmp_scaler, LinearRamp(0,100,0,1),r1=0.0, r2=0.0, r3=0.50, r4=0.70, k=250.0,
            atom_1_res_index=i, atom_1_name=a_i, atom_2_res_index=j, atom_2_name=a_j))

            hydroph_restraints.append(s.restraints.create_restraint_group(local_contact,1))
            #print 'hydroph:',i,j
    all_rest = len(hydroph_restraints)
    active = int( ContactsPerHydroph * len(group_1) )
    print(active,len(group_1),all_rest)
    s.restraints.add_selectively_active_collection(hydroph_restraints, active)


def create_HydrogenBond(s,HBPerResidue=0.10,scaler=None,group_1=np.array([]),group_2=np.array([]),CO=True):
    all_donor = ['N']
    all_acceptor = ['O']
    acceptors = {"SER":['OG'],
             "THR":['OG1'],
             "CYS":['SG'],
             "TYR":['OH'],
             "ASN":['OD1','ND2'],
             "GLN":['OE1','NE2'],
             "ASP":['OD1','OD2'],
             "GLU":['OE1','OE2'],
             "ARG":['NE','NH1','NH2'],
             "HIE":['ND1','NE2']}
    donors = {"SER":['OG'],
             "THR":['OG1'],
             "CYS":['SG'],
             "TYR":['OH'],
             "ASN":['ND2'],
             "GLN":['NE2'],
             "LYS":['NZ'],
             "ARG":['NE','NH1','NH2'],
             "HIE":['NE2'],
             "LYS":['NZ']}

    #Groups should be 1 centered
    n_res = s.residue_numbers[-1]
    group_1 = group_1 if group_1.size else np.array(list(range(n_res)))+1
    group_2 = group_2 if group_2.size else np.array(list(range(n_res)))+1
    scaler = scaler if scaler else s.restraints.create_scaler('nonlinear', alpha_min=0.4, alpha_max=1.0, factor=4.0)

    #Get a list of names and residue numbers, if just use names might skip some residues that are two
    #times in a row
    #make list 1 centered
    sequence = [(i,j) for i,j in zip(s.residue_numbers,s.residue_names)]
    sequence = sorted(set(sequence))
    sequence = dict(sequence)

    pairs = []
    HB_restraints = []
    for i in group_1:
        for j in group_2:

            # don't put the same pair in more than once
            if ( (i,j) in pairs ) or ( (j,i) in pairs ):
                continue

            if ( i ==j ):
                continue
            #Alpha helices are 1-4, exclude them from this since they have the ss heuristic
            if (abs(i-j)< 5):
                continue
            if (abs(i-j)> 10):
                continue
            pairs.append( (i,j) )

            if CO:
                tmp_scaler = scaler(abs(i-j),'Hbond')
            else:
                tmp_scaler = scaler
           
            local_contact = []
            #i as donor
            atoms_i = []
            atoms_i.extend(all_donor)
            try:  #some residues do not have donors
                atoms_i.extend(donors[sequence[i]])
            except:
                pass
            #j as acceptor
            atoms_j = []
            atoms_j.extend(all_acceptor)
            try:  #some residues do not have donors
                atoms_j.extend(acceptors[sequence[j]])
            except:
                pass

            #print "DONORS: ",i,sequence[i],atoms_i
            #print "ACCEPTORS: ",j,sequence[j],atoms_j

            for a_i in atoms_i:
                for a_j in atoms_j:
                    local_contact.append(s.restraints.create_restraint('distance', tmp_scaler, LinearRamp(0,100,0,1),r1=0.0, r2=0.0, r3=0.50, r4=0.65, k=250.0,
            atom_1_res_index=i, atom_1_name=a_i, atom_2_res_index=j, atom_2_name=a_j))

            #j as donor
            atoms_j = []
            atoms_j.extend(all_donor)
            try:  #some residues do not have donors
                atoms_i.extend(donors[sequence[j]])
            except:
                pass
            #i as acceptor
            atoms_i = []
            atoms_i.extend(all_acceptor)
            try:  #some residues do not have donors
                atoms_j.extend(acceptors[sequence[j]])
            except:
                pass
            for a_i in atoms_i:
                for a_j in atoms_j:
                    local_contact.append(s.restraints.create_restraint('distance', tmp_scaler, LinearRamp(0,100,0,1),r1=0.0, r2=0.0, r3=0.50, r4=0.65, k=250.0,
            atom_1_res_index=i, atom_1_name=a_i, atom_2_res_index=j, atom_2_name=a_j))

            HB_restraints.append(s.restraints.create_restraint_group(local_contact,1))
    all_rest = len(HB_restraints)
    active = int( HBPerResidue * len(group_1) )
    print(active,len(group_1),all_rest)
    s.restraints.add_selectively_active_collection(HB_restraints, active)

def gen_state(s, index):
    pos = s._coordinates
    pos = pos - np.mean(pos, axis=0)
    vel = np.zeros_like(pos)
    alpha = index / (N_REPLICAS - 1.0)
    s._box_vectors=np.array([0.,0.,0.])
    energy = 0
    return system.SystemState(pos, vel, alpha, energy,s._box_vectors)


def get_dist_restraints(filename, s, scaler):
    dists = []
    rest_group = []
    lines = open(filename).read().splitlines()
    lines = [line.strip() for line in lines]
    for line in lines:
        if not line:
            dists.append(s.restraints.create_restraint_group(rest_group, 1))
            rest_group = []
        else:
            cols = line.split()
            i = int(cols[0])
            name_i = cols[1]
            j = int(cols[2])
            name_j = cols[3]
            dist = float(cols[4]) / 10.

            rest = s.restraints.create_restraint('distance', scaler,LinearRamp(0,100,0,1),
                                                 r1=0.0, r2=0.0, r3=dist, r4=dist+0.2, k=250,
                                                 atom_1_res_index=i, atom_2_res_index=j,
                                                 atom_1_name=name_i, atom_2_name=name_j)
            rest_group.append(rest)
    return dists

def get_knob_restraints(filename, s, scaler):
    dists = []
    lines = open(filename).read().splitlines()
    lines = [line.strip() for line in lines]
    high = 0
    total = 0
    for line in lines:
            cols = line.split()
            i = int(cols[0])
            name_i = cols[2]
            j = int(cols[3])
            name_j = cols[5]
            dist = float(cols[6]) / 10.
            sd = float(cols[7]) / 10.
            accuracy = float(cols[8])

            if name_i not in ["CA","CB"]:
                continue
            if name_j not in ["CA","CB"]:
                continue
            if abs(i-j) < 6:
                continue
            if accuracy > 0.5:
                high = high + 1
                total = total + 1

            rest = s.restraints.create_restraint('distance', scaler,LinearRamp(0,100,0,1),
                                                 r1=max(0,dist-sd-0.2), r2=max(0,dist-sd), r3=dist+sd, r4=dist+sd+0.2, k=250,
                                                 atom_1_res_index=i, atom_2_res_index=j,
                                                 atom_1_name=name_i, atom_2_name=name_j)
            #collections are parallel; groups are in one gpu thread
            dists.append(s.restraints.create_restraint_group([rest], 1))
    return dists,float(high)/float(total)*0.95

def setup_system():
    # load the sequence
    sequence = parse.get_sequence_from_AA1(filename='sequence.dat')
    n_res = len(sequence.split())

    # build the system
    p = system.ProteinMoleculeFromSequence(sequence)
    b = system.SystemBuilder(forcefield="ff14sbside")
    s = b.build_system_from_molecules([p])
    s.temperature_scaler = system.GeometricTemperatureScaler(0, 0.6, 300., 450.)

    #
    # Secondary Structure
    #
    ss_scaler = s.restraints.create_scaler('constant')
    ss_rests = parse.get_secondary_structure_restraints(filename='ss.dat', system=s,ramp=LinearRamp(0,100,0,1), scaler=ss_scaler,
            torsion_force_constant=2.5, distance_force_constant=2.5)
    n_ss_keep = int(len(ss_rests) * 1.0) #We enforce 100% of restrains 
    s.restraints.add_selectively_active_collection(ss_rests, n_ss_keep)


    #
    # Distance Restraints
    #
    dist_scaler = s.restraints.create_scaler('nonlinear', alpha_min=0.4, alpha_max=1.0, factor=4.0)

    # High reliability
    #
    #
    old_protocol = s.restraints.create_scaler('nonlinear', alpha_min=0.40, alpha_max=1.00, factor=4.0)
    #
    # Heuristic Restraints
    #
    subset= np.array(list(range(n_res))) + 1

        #
        # Hydrophobic
        #
    create_hydrophobes(s,ContactsPerHydroph=1.2,scaler=old_protocol,group_1=subset,CO=False)
        #
        # Strand Pairing
        #
    sse,active = make_ss_groups(subset=subset)
    try:
        generate_strand_pairs(s,sse,float(active),subset=subset,scaler=old_protocol,CO=False)
    except:
        print("Not using Strand Pairing Heuristic")
        pass

        #
        # Evolutionary restraints
        #
    try:
        create_Evolution(s,accuracy=0.7,scaler=dist_scaler,fname='evolution_contacts.dat')
    except:
        print("Not using Evolutionary restraints")
        pass

    #
    # Distance Restraints
    #
        #
        # Knob restraints
        #
    try:
        knobs,knob_accuracy = get_knob_restraints('Knob.data',s,scaler=dist_scaler)
        n_knobs = int(len(knobs) * knob_accuracy) 
        s.restraints.add_selectively_active_collection(knobs,n_knobs)
    except:
        print("Not using Knob-Socket predictions")
        pass



    # setup mcmc at startup
    movers = []
    n_atoms = s.n_atoms
    for i in range(1, n_res + 1):
        n = s.index_of_atom(i, 'N') - 1
        ca = s.index_of_atom(i, 'CA') - 1
        c = s.index_of_atom(i, 'C') - 1

        mover = mc.DoubleTorsionMover(n, ca, list(range(ca, n_atoms)),
                                      ca, c, list(range(c, n_atoms)))

        movers.append((mover, 1))

    sched = mc.MonteCarloScheduler(movers, n_res * 60)

    # create the options
    options = system.RunOptions()
    options.implicit_solvent_model = 'gbNeck2'
    options.use_big_timestep = False
    options.use_bigger_timestep = True
    options.cutoff = 1.8

    options.use_amap = False
    options.amap_alpha_bias = 1.0
    options.amap_beta_bias = 1.0
    options.timesteps = 11111
    options.minimize_steps = 20000
   # for i in range(30):
   #     print("Heads up! using MC minimizer!")
#    options.min_mc = sched

    # create a store
    store = vault.DataStore(s.n_atoms, N_REPLICAS, s.get_pdb_writer(), block_size=BLOCK_SIZE)
    store.initialize(mode='w')
    store.save_system(s)
    store.save_run_options(options)

    # create and store the remd_runner
    l = ladder.NearestNeighborLadder(n_trials=100)
    policy = adaptor.AdaptationPolicy(2.0, 50, 50)
    a = adaptor.EqualAcceptanceAdaptor(n_replicas=N_REPLICAS, adaptation_policy=policy)

    remd_runner = master_runner.MasterReplicaExchangeRunner(N_REPLICAS, max_steps=N_STEPS, ladder=l, adaptor=a)
    store.save_remd_runner(remd_runner)

    # create and store the communicator
    c = comm.MPICommunicator(s.n_atoms, N_REPLICAS)
    store.save_communicator(c)

    # create and save the initial states
    states = [gen_state(s, i) for i in range(N_REPLICAS)]
    store.save_states(states, 0)

    # save data_store
    store.save_data_store()

    return s.n_atoms


setup_system()

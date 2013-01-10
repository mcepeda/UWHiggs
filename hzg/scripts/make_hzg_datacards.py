#!/usr/bin/env python

from UWHiggs.hzg.datacard.directory_prep import directory_prep
from UWHiggs.hzg.datacard.metadata_association import metadata_association
from UWHiggs.hzg.datacard.categories_map import categories_map

import os

import ROOT
from ROOT import TFile,TH1F,TTree,RooWorkspace,gDirectory,RooDataSet,RooFit,\
     RooArgSet, RooFormulaVar, RooArgList, RooConstVar

def initialize_workspace(ws):
    #make the category definition
    cattype   = os.environ['hzgcategorytype']
    catstring = categories_map[cattype]['catstring']
    catname   = categories_map[cattype]['leafname']
    ws.factory(catstring)

    #fitting variables (init with no implicit range)
    ws.factory("Mzg[0]")
    
    #(possible) weights
    ws.factory('procWeight[0]')
    ws.factory('puWeight[0]')
    ws.factory('weight[0]')

    ws.defineSet('vars_with_weights','Mzg,procWeight,puWeight,%s'%catname)
    ws.defineSet('vars_with_weights_final',
                 'Mzg,procWeight,puWeight,weight,%s'%catname)
    ws.defineSet('vars','Mzg,%s'%catname)

def make_weighted_dataset(subproc,ws,tree,mc_events):
    data = RooDataSet('%s_shape_data'%subproc,
                      'M_{ll#gamma} Shape Data for %s'%subproc,
                      tree,
                      ws.set('vars_with_weights')                      
                      )

    mc_yield_var = RooConstVar('temp','temp',mc_events)
    weighter = RooFormulaVar('weight','weight','@0*@1/@2',
                             RooArgList( ws.var('procWeight'),
                                         ws.var('puWeight'),
                                         mc_yield_var )
                             )
    data.addColumn(weighter)

    data_proper_weight = RooDataSet('%s_shape_data'%subproc,
                                    'M_{ll#gamma} Shape Data for %s'%subproc,
                                    data,
                                    ws.set('vars_with_weights_final'),
                                    '','weight')
    return data_proper_weight

def extract_higgs_data_in_categories(subproc,input_file,ws):
    pwd = gDirectory.GetPath()    
    fin = TFile.Open(input_file,'read')
    gDirectory.cd(pwd)

    tree = fin.Get('selected_zg')    

    mc_tot_events = float(fin.Get('eventCount').GetBinContent(1))

    data_proper_weight = make_weighted_dataset(subproc,ws,tree,mc_tot_events)
    
    fin.Close()    
    
    getattr(ws,'import')(data_proper_weight)

def extract_bkg_data_in_categories(subproc,input_file,ws):
    pwd = gDirectory.GetPath()    
    fin = TFile.Open(input_file,'read')
    gDirectory.cd(pwd)

    tree = fin.Get('selected_zg')    

    mc_tot_events = float(fin.Get('eventCount').GetBinContent(1))

    data_proper_weight = make_weighted_dataset(subproc,ws,tree,mc_tot_events)
    
    fin.Close()    
    
    data_in_ws = ws.data('background_shape_data')
    if not not data_in_ws:
        getattr(ws,'import')(data_proper_weight)
        data_in_ws.append(data_proper_weight)        
    else:
        getattr(ws,'import')(data_proper_weight,
                             RooFit.Rename('background_shape_data'))
        getattr(ws,'import')(data_proper_weight)
                             

def extract_data_in_categories(channel,input_file,ws):
    pwd = gDirectory.GetPath()    
    fin = TFile.Open(input_file,'read')
    gDirectory.cd(pwd)

    tree = fin.Get('selected_zg')

    data = RooDataSet('%s_data'%channel,
                      'real data %s channel'%channel,
                      tree,
                      ws.set('vars') )    


    data_in_ws = ws.data('%s_data'%channel)
    if not not data_in_ws:
        data_in_ws.append(data)
    else:
        getattr(ws,'import')(data)

def create_master_workspaces(meta_data):
    pwd = gDirectory.GetPath()
    ws_list = {}
    for (sample,chanlist) in meta_data.getAssociation().iteritems():
        ws_list[sample] = {}
        for (channel,proclist) in chanlist.iteritems():            
            #make the workspace we're going to use
            this_ws = RooWorkspace('%s-%s'%(channel,sample))
            initialize_workspace(this_ws)
            ws_list[sample][channel] = {}
            for (process,subproclist) in proclist.iteritems():
                print sample, channel, process
                if 'HToZG'in process:
                    ws_list[sample][channel][process] = []
                for (subproc,info) in subproclist.iteritems():
                    print '\tprocessing: %s'%subproc
                    input_file = info['input_file']
                    if input_file == '':
                        print '\t no input file found! Skipping!'
                        continue
                    if 'data' not in process:                        
                        print '\t mc input = %s'%input_file.split('/')[-1]
                        if 'HToZG' in subproc:
                            ws_list[sample][channel][process].append(subproc)
                            extract_higgs_data_in_categories(subproc,
                                                             input_file,
                                                             this_ws)
                        else:
                            extract_bkg_data_in_categories(subproc,
                                                           input_file,
                                                           this_ws)
                    else:
                        print '\t data input = %s'%input_file.split('/')[-1]
                        extract_data_in_categories(channel,input_file,this_ws)
            #end loop over processes and data
            fout_name = '%s_%s_master_workspace.root'%(channel,sample)
            fout = TFile.Open(fout_name,'recreate')
            fout.cd()
            this_ws.Write()
            fout.Close()
            gDirectory.cd(pwd)
            ws_list[sample][channel]['filename'] = fout_name        
    return ws_list

def build_category_workspaces(ws_list,metadata):
    assoc = metadata.getAssociation()
    for (sample,chanlist) in ws_list.iteritems():
        for (channel,chaninfo) in chanlist.iteritems():
            processes = chaninfo.keys()
            processes.remove('filename')
            for process in processes:
                subproc_list = chaninfo[process]
                print sample, channel, process, subproc_list

if __name__ == '__main__':
    analysis_root = os.environ['hzganalysisroot']
    analysis_name = os.environ['hzganalysisname']
    
    dp = directory_prep(analysis_root,analysis_name)
    dp.build_inputs()

    mda = metadata_association(dp)
    
    master_ws_list = create_master_workspaces(mda)

    build_category_workspaces(master_ws_list,mda)
                
    


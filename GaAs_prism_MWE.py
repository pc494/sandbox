import pyprismatic as pr

meta_params = {}
meta_params['save4DOutput'] = True
meta_params['save3DOutput'] = False
meta_params['scanWindowXMin'] = 0.49
meta_params['scanWindowXMax'] = 0.51
meta_params['scanWindowYMin'] = 0.49
meta_params['scanWindowYMax'] = 0.51
meta_params['E0']=100

meta = pr.Metadata(filenameAtoms="PP_input.XYZ",filenameOutput="PP_output_GaAs",**meta_params)
meta.go()
output = pr.fileio.readMRC(meta.filenameOutput+'X0_Y0_FP1.mrc')
print(output)

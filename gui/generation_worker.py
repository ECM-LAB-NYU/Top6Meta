# =============================================================================
# Top-6-Class MetaStudio (Top6Meta)
# Version: 1.0
#
# A Python HPC framework for modeling architected materials and metastructures.
#
# Authors:
#   - Agyapal Singh [1]
#   - Georgios Mermigkis [2]
#   - Panagiotis Hadjidoukas [2]
#   - Nikolaos Karathanasopoulos [1]
#
# Affiliations:
#   [1] New York University, Department of Engineering,
#       Abu Dhabi, United Arab Emirates
#   [2] Computer Engineering and Informatics Department,
#       University of Patras, Greece
#
# © 2026 The Authors
#
# License: MIT License
# =============================================================================
import sys
import os
import pickle
import traceback

# parent directory to path to import tpms_core
sys.path.append('../python')
import tpms_core

def stop_callback():
    # This callback can later be used to check if a user wants to stop the generation
    # process earlier
    return False

def main():
    if len(sys.argv) != 3:
        print("Usage: python generation_worker.py <params_file> <temp_dir>", file=sys.stderr)
        sys.exit(1)
    
    params_file = sys.argv[1]
    temp_dir = sys.argv[2]
    
    try:
        try:
            import os
            os.nice(-5)  # Higher priority on Unix systems
        except:
            pass
        
        # load parameters
        with open(params_file, 'rb') as f:
            option, params = pickle.load(f)
        
        print(f"Starting {option} generation...", file=sys.stderr)
        
        # run the computation based on option
        if option == "TPMS":
            if params.get('IPC') == "IPC_Y":
                Freinf, Vreinf, Fcompl, Vcompl, Final_Vol_Frac, Final_Surface = tpms_core.generate_tpms(
                    **params, stop_callback=stop_callback
                )
                result = {
                    "F_reinf": Freinf,
                    "V_reinf": Vreinf,
                    "F_compl": Fcompl,
                    "V_compl": Vcompl,
                    "Final_Vol_Frac": Final_Vol_Frac,
                    "Final_Surface": Final_Surface
                }
            else:
                F, V, Final_Vol_Frac, Final_Surface = tpms_core.generate_tpms(
                    **params, stop_callback=stop_callback
                )
                result = {
                    "F": F,
                    "V": V,
                    "Final_Vol_Frac": Final_Vol_Frac,
                    "Final_Surface": Final_Surface
                }
            
        elif option == "Spinodal":
            if params.get('IPC') == "IPC_Y":
                Freinf, Vreinf, Fcompl, Vcompl, Final_Vol_Frac, Final_Surface = tpms_core.generate_tpms(**params)
                result = {
                    "F_reinf": Freinf,
                    "V_reinf": Vreinf,
                    "F_compl": Fcompl,
                    "V_compl": Vcompl,
                    "Final_Vol_Frac": Final_Vol_Frac,
                    "Final_Surface": Final_Surface
                }
            else:
                F, V, Final_Vol_Frac, Final_Surface = tpms_core.generate_tpms(**params)
                result = {
                    "F": F,
                    "V": V,
                    "Final_Vol_Frac": Final_Vol_Frac,
                    "Final_Surface": Final_Surface
                }
            
        elif option == "Strut":
            if params.get('IPC') == "IPC_Y":
                Freinf, Vreinf, Fcompl, Vcompl, Final_Vol_Frac, Final_Surface = tpms_core.generate_strut(**params)
                result = {
                    "F_reinf": Freinf,
                    "V_reinf": Vreinf,
                    "F_compl": Fcompl,
                    "V_compl": Vcompl,
                    "Final_Vol_Frac": Final_Vol_Frac,
                    "Final_Surface": Final_Surface
                }
            else:
                F, V, Final_Vol_Frac, Final_Surface = tpms_core.generate_strut(**params)
                result = {
                    "F": F,
                    "V": V,
                    "Final_Vol_Frac": Final_Vol_Frac,
                    "Final_Surface": Final_Surface
                }
            
        elif option == "Hybrid":
            if params.get('IPC') == "IPC_Y":
                Freinf, Vreinf, Fcompl, Vcompl, Final_Vol_Frac, Final_Surface = tpms_core.generate_hybrid(**params)
                result = {
                    "F_reinf": Freinf,
                    "V_reinf": Vreinf,
                    "F_compl": Fcompl,
                    "V_compl": Vcompl,
                    "Final_Vol_Frac": Final_Vol_Frac,
                    "Final_Surface": Final_Surface
                }
            else:
                F, V, Final_Vol_Frac, Final_Surface = tpms_core.generate_hybrid(**params)
                result = {
                    "F": F,
                    "V": V,
                    "Final_Vol_Frac": Final_Vol_Frac,
                    "Final_Surface": Final_Surface
                }
            
        elif option == "Layered":
            F0, V0, F1, V1, Final_Vol_Frac, Final_Surface = tpms_core.generate_layered(**params)
            result = {
                "F0": F0,
                "V0": V0,
                "F1": F1,
                "V1": V1,
                "Final_Vol_Frac": Final_Vol_Frac,
                "Final_Surface": Final_Surface
            }
            
        else:
            raise ValueError(f"Invalid option for generation: {option}")
        
        # Save result
        result_file = os.path.join(temp_dir, "result.pkl")
        with open(result_file, 'wb') as f:
            pickle.dump(result, f)
        
        print(f"{option} generation completed successfully!", file=sys.stderr)
        
    except Exception as e:
        # Save error
        error_file = os.path.join(temp_dir, "error.txt")
        with open(error_file, 'w') as f:
            f.write(f"Error in {option} generation: {str(e)}\n")
            f.write(traceback.format_exc())
        
        print(f"Error in {option} generation: {str(e)}", file=sys.stderr)
        sys.exit(1)
    
    finally:
        # clean up params file
        try:
            os.remove(params_file)
        except:
            pass

if __name__ == "__main__":
    main()
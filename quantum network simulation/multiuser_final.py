

"""
Improved QKD Network Implementation with Wavelength Division Multiplexing
"""

import SingleLink as sl
import numpy as np
from typing import Tuple, List, Dict
from scipy.optimize import differential_evolution
import logging
import random
import csv
import os



class QKDNetwork:
    def __init__(self, num_users: int):
        self.num_users = num_users
        
 
        self.wa_subnet1 = np.array([
            [5, -5, -3, 5, -5, -3],
            [4, 3, -4, 4, 3, -4],
            [2, 1, 1, -2, -1, -1]])
        
        self.wa_subnet2 = np.array(
           [[6,-6,-4,6,-6,-4],
            [5,4,-5, 5,4,-5],
            [3,2,1,-3,-2,-1]])
        

        self.wa_fullmesh=np.array(
            [[1,-1,-2,-3,-4,-5],
             [2,6,-6,-7,-8,-9],
             [3,7,10,-10,-11,-12],
             [4,8,11,13,-13,-14],
             [5,9,12,14,15,-15]])
        
        self.wa_add_channel1=np.array(
            [[1,-1,-2,-3,-4,-5],
             [2,6,-6,-7,-8,-9],
             [3,7,10,-10,-11,-12],
             [4,8,11,13,-13,-14],
             [5,9,12,14,15,-15],
             [16,-16,0,0,0,0]])
        
        self.wa_add_channel2=np.array(
            [[1,-1,-2,-3,-4,-5],
             [2,6,-6,-7,-8,-9],
             [3,7,10,-10,-11,-12],
             [4,8,11,13,-13,-14],
             [5,9,12,14,15,-15],
             [16,-16,0,0,0,0],
             [17,-17,0,0,0,0]])
        
        #AB AC are not connected
        self.wa_partial_mesh1=np.array(
            [[0,0,0,-1,-2,-3],
             [0,4,-4,-5,-6,-7],     
             [1,5,8,-8,-9,-10],
             [2,6,9,11,-11,-12],
             [3,7,10,12,13,-13]])
        
        #AD AF AG BC are not connected
        self.wa_partial_mesh2=np.array(
            [[1,-1,-2,-3,-4,-5],
             [2,3,6,-6,-7,-8],
             [0,4,7,9,-9,-10],
             [0,5,8,10,11,-11]])
        
        
       
    
    #check if the wavelength_allocation is valid
    def validate_wavelength_allocation(self,wavelength_allocation: np.ndarray)-> np.ndarray: 
        sender=[]
        receiver=[]
        num_cols = wavelength_allocation.shape[1]  
        sender, receiver = np.triu_indices(num_cols, k=1)  # Get upper triangle indices
        
        for i in range(len(sender)):
            chan1 = wavelength_allocation[:, sender[i]].tolist()  
            chan2 = wavelength_allocation[:, receiver[i]].tolist()
        
        #check if all users have channels to use
        if not any(-item in chan2 for item in chan1):
            print('warning: some links are not connected')

        return wavelength_allocation
    
    
    '''
    return the name of wavelength allocation at the output stage
    '''
    def find_wavelength_allocation(self, wavelength_allocation: np.ndarray) -> str:
        wavelength_allocation_map = {
        'subnet1': self.wa_subnet1,
        'subnet2': self.wa_subnet2,
        'fullmesh': self.wa_fullmesh,
        'additional_chan1': self.wa_add_channel1,
        'additional_chan2': self.wa_add_channel2,
        'partial_mesh1': self.wa_partial_mesh1,
        'partial_mesh2': self.wa_partial_mesh2
        }

        for name, value in wavelength_allocation_map.items():
            if np.array_equal(wavelength_allocation, value):
                return name
    
    '''
    generate random distance from user to QNSP
    '''
    def generate_fiber_length(self, wa:int, num_users: int)->np.ndarray:
        if wa==0:
            fiber_length=[random.randint(1, 30) for _ in range(num_users)]
            #fiber_length=[15, 3, 29, 28, 17, 29]
            self.fb_variable=fiber_length
        if wa>0:
            fiber_length=self.fb_variable
        #fiber_length=[26, 30, 2, 27, 22, 8]
        #print(fiber_length)
        return fiber_length
        
        
                
                

    '''
    (1)assume no counts from other users and detector efficiency is 70%, total pair 
    generation rate=10^6
    (2)input: users: order of this user in the sequence E.g. Alice->0, Chloe->2
    (3)output: ideal_single_link: counts at the link correspond to the user assume counts
    from other user=0
    pgr_per_channel: pair generation rate per cahnnel corresponds to wavelength allocation
    '''

    def ideal_single_counts(self, wavelength_allocation: np.ndarray, user: int,\
                          fiber_length: np.ndarray, pair_generation_rate:int, bs=int) -> Tuple[int, int]:
 
        ideal_single_link = np.zeros(self.num_users, dtype=int)
        pgr_per_channel = (pair_generation_rate) / np.max(wavelength_allocation)
        
        singlecounts=sl.Calc_QKD_Parameters_For_Single_Link\
        (PairGenerationRate_per_second=pgr_per_channel, FibreLength_in_km = \
        [fiber_length[user],0],SingleCountsFromOtherUsers_per_second=[0,0,0,0], 
        OtherLoss_in_dB = [6+(bs-1)*3,6+(bs-1)*3])
        
        ideal_single_link = singlecounts['Singles rate for Alice'].HA / 0.7
        return ideal_single_link, pgr_per_channel
    
    
    
    '''
    (1)This function calculates the number of useless channels (noise) and 
    effective channels for the specific node pairs
    (2)input: row and col: sender and receiver nodes
    (3)output: channels_at_nodes: availble channels for both sender and receiver
               effective_channel: channels can be used for effective communication
               num_eff_chan: number of effective channel at the sender node
               num_noise_chan: number of noise channel at the sender node
               'full' channel: this pair of wavelength is used only once
               'half' channel: this pair of wavelength is used twice (only half keys are valid)
    '''

    def num_channel_from_other_users(self, wavelength_allocation: np.ndarray,
                                   row: int, col: int) -> Tuple[List, List, int, int]:

        channels_send = []
        channels_receive= []
        effective_channel = []
        count_dict = {}
        
        
        num_full_eff = 0
        num_half_eff = 0
        num_full_noise = 0
        num_half_noise = 0

         #first half of the list is row and second half is col
        for x in range(len(wavelength_allocation)):
            for y in range(len(wavelength_allocation[0])):
                if y == row and wavelength_allocation[x][y]!=0:
                    channels_send.append(wavelength_allocation[x][y])
                    
        for x in range(len(wavelength_allocation)):
            for y in range(len(wavelength_allocation[0])):
                if y == col and wavelength_allocation[x][y]!=0:
                    channels_receive.append(wavelength_allocation[x][y])
            
        
        for num in channels_send:
            if -num in channels_receive :
                effective_channel.append(abs(num))
        #print(effective_channel)
        #record the number of wavelength occured at the wavelength assignment    
        for x in wavelength_allocation:
            for y in x:
                count_dict[y]=count_dict.get(y, 0) + 1       
        for element in effective_channel:
            if count_dict.get(element, 0)==1:
                num_full_eff += 1
            if count_dict.get(element, 0)==2:
                num_half_eff += 1
        for element in channels_send:
            if count_dict.get(-element, 0)==1:
                num_full_noise += 1
            if count_dict.get(-element, 0)==2:
                num_half_noise += count_dict.get(-element, 0)
        num_full_noise=num_full_noise-num_full_eff
        num_half_noise=num_half_noise-num_half_eff
      
        return num_full_eff, num_half_eff, num_full_noise, num_half_noise
    
    
    
    '''
    calculate the total noise counts from other users given with sender and receiver
    '''
    def pgr_noise_counts(self, wavelength_allocation: np.ndarray, fiber_length: np.ndarray, 
                          pair_generation_rate:int, row: int, col: int) -> Tuple[int, int]:
     
        num_full_eff, num_half_eff, num_full_noise, num_half_noise = self.num_channel_from_other_users(
            wavelength_allocation, row, col
            )
        single_counts_full, pgr = self.ideal_single_counts(
            wavelength_allocation, row, fiber_length, pair_generation_rate, 1
            )
        single_counts_half, _ = self.ideal_single_counts(
            wavelength_allocation, row, fiber_length, pair_generation_rate, 2
            )
        link_pgr = pgr * (num_full_eff+num_half_eff)
        node_noise = single_counts_full * num_full_noise + single_counts_half * num_half_noise
        
        if num_half_eff==0:
            bs_loss=0
        else:
            if num_full_eff==0:
                bs_loss=3

        return link_pgr, node_noise, bs_loss


    '''
    calculate the new secure key rate of each link base on the counts from other user
    and the number of effective channel 
    '''
    def multi_secure_key_rate(self, wavelength_allocation: np.ndarray,fiber_length:np.ndarray,
                 num_users: int, coincidence_matrix: np.ndarray, pair_generation_rate: int) -> Tuple[List, List]:
        
        secure_key_rate_matrix = np.zeros((self.num_users, self.num_users), dtype=int)
        QBER_matrix=np.zeros((self.num_users, self.num_users), dtype=float)
        for row in range(self.num_users):
            for col in range(self.num_users):
                if row > col:
                        
                    new_pgr, node_noise1, bs_loss1=self.pgr_noise_counts(wavelength_allocation, 
                                                fiber_length, pair_generation_rate, row, col)
                    _, node_noise2, bs_loss2=self.pgr_noise_counts(wavelength_allocation, 
                                                fiber_length, pair_generation_rate, col, row)
                    new_single_link = sl.Calc_QKD_Parameters_For_Single_Link(
                        PairGenerationRate_per_second=new_pgr,
                        CoincidenceWindow_in_ps=coincidence_matrix[row][col],
                        FibreLength_in_km=[fiber_length[row], fiber_length[col]],
                        OtherLoss_in_dB = [6+bs_loss1,6+bs_loss2],
                        SingleCountsFromOtherUsers_per_second=[
                            node_noise1, node_noise1,
                            node_noise2, node_noise2
                            ]
                    )
                
                    secure_key_rate_matrix[row][col] = new_single_link['Total secure key rate']
                    secure_key_rate_matrix[col][row] = new_single_link['Total secure key rate']
                           
                    QBER_matrix[row][col]=new_single_link['QBER between Alice and Bob']
                    QBER_matrix[col][row]=new_single_link['QBER between Alice and Bob']
        sum_of_skr = np.sum(secure_key_rate_matrix)/2
        min_of_skr = min(self.matrix_to_list(secure_key_rate_matrix))
        min_parameter=min_of_skr*3000+sum_of_skr
        #if min_of_skr<1:
            #return sum_of_skr-9999, min_of_skr-99, min_parameter-9999, secure_key_rate_matrix, secure_key_rate_matrix, QBER_matrix
        return sum_of_skr, min_of_skr, min_parameter, secure_key_rate_matrix, QBER_matrix
    
    
    def list_to_matrix(self, params: np.ndarray) -> np.ndarray:
        matrix = np.zeros((self.num_users, self.num_users))
        idx = 0
        for i in range(self.num_users):
            for j in range(i + 1, self.num_users):
                matrix[i,j] = matrix[j,i] = params[idx]
                idx += 1
        return matrix
    
    def matrix_to_list(self, params: np.ndarray) -> np.ndarray:
        plotlist = []
        for row in range(len(params)):
            for col in range(len(params[0])):
                if row<col:
                    plotlist.append(params[row][col])
        return plotlist
        
    '''
    (1)pgr_flag=0: optimized pgr, pgr_flag=1: fixed pgr
    (2)opt_mode=0: total skr priority, opt=1: minimum skr priority
    (3)use differential evolution to find the best combination of PGR and CW, use the best parameters 
    as input to revaluate the performance of the network
    '''
    def optimize_network_parameters(self, wavelength_allocation: np.ndarray, num_users:int, wa:int,
                                    pgr_flag:int, opt_mode:int) -> Dict:
        

        num_channel = (self.num_users * (self.num_users - 1)) // 2 
        fiber_length=self.generate_fiber_length(wa, num_users)
        print('fiber length：',fiber_length)
        #set boundary for coincidence window and pair generation rate
        bounds = [(40, 800)] * num_channel + [(4e2, 1e4)]
        
        if pgr_flag==0:
            wrapped_objective = lambda x: -self.multi_secure_key_rate(wavelength_allocation,fiber_length, num_users,
                                                                  self.list_to_matrix(x[:-1]), x[-1]*1000 )[opt_mode]
        if pgr_flag==1:
            wrapped_objective = lambda x: -self.multi_secure_key_rate(wavelength_allocation,fiber_length, num_users,
                                                                      self.list_to_matrix(x[:-1]), 1e6)[opt_mode]
        
        base_popsize = 35
        max_attempts = 5
        attempt = 0
        successful = False
        while not successful and attempt < max_attempts:
            popsize = base_popsize + attempt * 5
            print(f"Attempt {attempt+1}: popsize={popsize},") 
            result = differential_evolution(
                wrapped_objective,
                bounds,
                strategy='rand1bin',
                popsize=popsize,
                maxiter=1000,
                mutation=(0.8,1.7),
                tol=0.008,
                recombination=0.9,
                polish=False,
                disp=True)
            
            if result.fun < 0 and result.nit > 1:
                successful = True
            else:
                attempt += 1
        
        if not result.success:
            logging.warning(f"Optimization did not converge: {result.message}")
        
        
        best_params = result.x
        n_cw = (self.num_users * (self.num_users - 1)) // 2
        best_cw_params = best_params[:n_cw]
        if pgr_flag==0:
            best_pgr = best_params[n_cw]*1000
        if pgr_flag==1:
            best_pgr = 1e6
        
        best_cw = np.zeros((self.num_users, self.num_users))
        idx = 0
        for i in range(self.num_users):
            for j in range(i + 1, self.num_users):
                best_cw[i,j] = best_cw[j,i] = int(best_cw_params[idx])
                idx += 1


        
        sum_skr,min_skr,min_parameter, skr_matrix,qber_matrix = self.multi_secure_key_rate(wavelength_allocation, fiber_length,
                     num_users, best_cw, best_pgr)
            
        average_qber=np.sum(qber_matrix)/(self.num_users*(self.num_users-1))   
        print(skr_matrix)
        if min_parameter<0:
            return None
        
                    
        wavelength_allocation_order=self.find_wavelength_allocation(wavelength_allocation)
        
       
        
        return {
            'coincidence_window': best_cw,
            'pair_generation_rate': best_pgr,
            'sum of secure key rate': sum_skr,
            'secure key rate matrix': skr_matrix,
            'average qber': average_qber,
            'QBER matrix': qber_matrix,
            'optimization_success': result.success,
            'optimization_message': result.message,
            'number_of_iterations': result.nit,
            'fiber_length': fiber_length,
            'wavelength_allocation': wavelength_allocation_order
        }
    
    
    '''
    create the CSV file if it's not existed in the folder, set the name of each column
    '''    
    def check_and_create_csv(self,filename='output.csv'):
        data = ["fiber_length", "total secure key rate","average QBER", "coincidence window",
                "pair generation rate", "number of iterations","optimization_success",
                "wavelength allocation order"]
        if not os.path.exists(filename):
            with open(filename, mode='w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(data)
    
    def write_csv(self, results):
        
        
        csv_store = [
            results['fiber_length'],  
            int(results['sum of secure key rate']),  
            round(float(results['average qber']), 4),  
            [int(param) for param in self.matrix_to_list(results['coincidence_window'])],
            round(results['pair_generation_rate'] / 1000) * 1000,  
            results['number_of_iterations'],  
            bool(results['optimization_success']),
            results['wavelength_allocation']
            ]
        with open("output4.csv", mode="a", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(csv_store)
        print("New data added to CSV file")            
        
        
'''
(1)add the name of wavelength allocation to be tested at wa_var
(2)change the last two imput parameters as flags at network.optimize_network_parameters for 
different optimization target and pair generation rate setting 
(3)fiber length combination would be the same among all the wavelength allocation scheme in wa_var
'''

def main():

    network = QKDNetwork(num_users=6)
    
  
    logging.basicConfig(level=logging.INFO)
    wa_var=[network.wa_add_channel2,network.wa_subnet1, network.wa_subnet2]
    
    network.check_and_create_csv('output4.csv')
    for m in range(len(wa_var)):
        i=m%len(wa_var)
        
        wavelength_allocation = network.validate_wavelength_allocation(wa_var[i])
        
        logging.info("Starting network optimization...")

        results = network.optimize_network_parameters(wavelength_allocation, network.num_users,i,0,0)
        network.write_csv(results)
        
        
        
        logging.info("\nOptimization Results:")
        logging.info(f"Best pair generation rate: {results['pair_generation_rate']:.2e}")
        logging.info(f"secure key rate matrix: \n {results['secure key rate matrix']}")
        logging.info(f"total secure key rate: {results['sum of secure key rate']:.2f}")
        logging.info(f"Average QBER: {results['average qber']:.4f}")
        logging.info(f"Number of iterations: {results['number_of_iterations']}")
       
        
    
if __name__ == "__main__":
    main()
    
    
    
    
    
    

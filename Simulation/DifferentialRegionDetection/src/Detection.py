import numpy as np
import cvxpy as cp
from scipy.sparse import csr_matrix
import logging


class EventDetection:
    def __init__(self, shorter_genes_ratio: list, samples_condition1, lambda2: list):
        self._shorter_genes_ratio = shorter_genes_ratio
        self._samples_condition1 = samples_condition1
        self._lambda2 = lambda2

    def changepoint_position(self):
        la2 = self._lambda2
        beta = EventDetection.fusedlasso(self._shorter_genes_ratio, la2)
        cp_position_index = EventDetection.search_cp(beta)
        return beta, cp_position_index

    @staticmethod
    def fusedlasso(data: list, lambda2):
        m = len(data)
        # Define array G
        row_a = np.array(range(m - 1))
        col_a = np.array(range(m - 1))
        data_a = np.array([-1] * (m - 1))
        a = csr_matrix((data_a, (row_a, col_a)), shape=(m - 1, m))

        row_b = np.array(range(m - 1))
        col_b = np.array(range(1, m))
        data_b = np.array([1] * (m - 1))
        b = csr_matrix((data_b, (row_b, col_b)), shape=(m - 1, m))

        G = (a + b).toarray()

        # create optimization variables.
        beta = cp.Variable(m)
        cost = cp.sum_squares(data - beta)
        obj = cp.Minimize(cost + lambda2 * cp.norm(G * beta, 1))
        # Form and solve problem.
        prob = cp.Problem(obj)
        prob.solve(solver='OSQP')
        return beta.value

    @staticmethod
    def search_cp(data: list):
        cp_index = []
        for i in range(len(data)-1):
            if abs(data[i]-data[i+1]) >= 0.01:
                cp_index.append(i+1)
        if len(cp_index) > 0:
            return cp_index
        else:
            return None

    @staticmethod
    def decide_close_cp(cp_position: list):
        n = len(cp_position)
        decide = []
        for i in range(n-1):
            if abs(cp_position[i]-cp_position[i+1]) <= 50:
                decide.append(1)
            else:
                decide.append(0)
        if decide.count(0) == 0:
            return True
        else:
            return False

    @staticmethod
    def get_ratio_and_dk(k_list, lam_list, j_list):
        if len(k_list) <= 5:
            logging.debug('data too small')
            return None
        k_list.reverse()
        lam_list.reverse()
        j_list.reverse()
        n = len(lam_list)

        lam_i_1 = [np.nan] + lam_list[:-1]

        li = [lam_i_1[i]-lam_list[i] for i in range(n)]

        ratio = []
        subb = []
        for i in range(n-1):
            if li[i] == np.nan:
                ratio.append(np.nan)
                subb.append(np.nan)
            else:
                get_max = max(li[i+1:])
                ratio.append(li[i]/get_max)
                subb.append(li[i]-get_max)
        ratio.append(np.nan)
        subb.append(np.nan)

        jk = [(j_list[-1]-j_list[i]) * (k_list[-1]-1) / (j_list[-1]-j_list[0])+1 for i in range(n)]

        dk = [np.nan] + [jk[i-1]-2*jk[i]+jk[i+1] for i in range(1, n-1)] + [np.nan]
        return ratio, subb, dk, li

    @staticmethod
    ######待定#######
    def determine_lambda2_and_cp_num(k_list, lam_list, ratio_list, sub_list, dk_list, thres):
        lambda2_final = None
        number_cp = None
        # determination
        max_ratio_thre = max(ratio_list[1:len(ratio_list) - 3])
        ratio_sub = list(map(lambda x, y: x * y, ratio_list, sub_list))
        max_ratio = max(ratio_sub[1:len(ratio_sub) - 2])
        max_index = ratio_sub.index(max_ratio)
        max_dk = max(dk_list[1:len(dk_list) - 1])
        max_dk_index = dk_list.index(max_dk)
        if max_ratio_thre >= 2 and max_dk > 5:
            logging.debug('dk is large')
            lambda2_final = lam_list[max_dk_index] + 1
            number_cp = k_list[max_dk_index] - 1
            return lambda2_final, number_cp
        if max_dk < 1.5:
            logging.debug('max dk is small, so change threshold to 1')
            thres = 1
        if max(dk_list[max_index + 1:len(dk_list) - 1]) > thres:
            logging.debug('k>kmax的dk有大于{}, and max dk is {}'.format(thres, max_dk))
            lambda2_final = lam_list[0] + 1
            number_cp = 0
            return lambda2_final, number_cp
        else:
            if dk_list[max_index] > thres:
                logging.debug('k=kmax, and dk is {}'.format(dk_list[max_index]))
                lambda2_final = lam_list[max_index] + 1
                number_cp = k_list[max_index] - 1
                return lambda2_final, number_cp
            else:
                if max_index == 1:
                    logging.debug('kmax=1并且dk<{}'.format(thres))
                    lambda2_final = lam_list[0] + 1
                    number_cp = 0
                    return lambda2_final, number_cp
                else:
                    for i in range(1, max_index):
                        iindex = max_index - i
                        if dk_list[iindex] > thres:
                            logging.debug('k在kmax上面 and dk is {}'.format(dk_list[iindex]))
                            lambda2_final = lam_list[iindex]+1
                            number_cp = k_list[iindex]-1
                            return lambda2_final, number_cp
                        else:
                            continue
                    if lambda2_final is None or number_cp is None:
                        return lam_list[0]+1, 0
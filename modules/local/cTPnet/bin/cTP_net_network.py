import torch.nn as nn
import torch.nn.functional as F
import torch

"""
Structure of the 12 MB-DNN for cTPnet.
"""
class Net12(nn.Module):
    def __init__(self):
        super(Net12, self).__init__()
        self.fc1 = nn.Linear(12611, 1000)
        self.fc3 = nn.Linear(1000, 128)
        self.fc40 = nn.Linear(128, 64)
        self.fc50 = nn.Linear(64, 1)
        self.fc41 = nn.Linear(128, 64)
        self.fc51 = nn.Linear(64, 1)
        self.fc42 = nn.Linear(128, 64)
        self.fc52 = nn.Linear(64, 1)
        self.fc43 = nn.Linear(128, 64)
        self.fc53 = nn.Linear(64, 1)
        self.fc44 = nn.Linear(128, 64)
        self.fc54 = nn.Linear(64, 1)
        self.fc45 = nn.Linear(128, 64)
        self.fc55 = nn.Linear(64, 1)
        self.fc46 = nn.Linear(128, 64)
        self.fc56 = nn.Linear(64, 1)
        self.fc47 = nn.Linear(128, 64)
        self.fc57 = nn.Linear(64, 1)
        self.fc48 = nn.Linear(128, 64)
        self.fc58 = nn.Linear(64, 1)
        self.fc49 = nn.Linear(128, 64)
        self.fc59 = nn.Linear(64, 1)
        self.fc410 = nn.Linear(128, 64)
        self.fc510 = nn.Linear(64, 1)
        self.fc411 = nn.Linear(128, 64)
        self.fc511 = nn.Linear(64, 1)
    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc3(x))
        o0 = self.fc50(F.relu(self.fc40(x)))
        o1 = self.fc51(F.relu(self.fc41(x)))
        o2 = self.fc52(F.relu(self.fc42(x)))
        o3 = self.fc53(F.relu(self.fc43(x)))
        o4 = self.fc54(F.relu(self.fc44(x)))
        o5 = self.fc55(F.relu(self.fc45(x)))
        o6 = self.fc56(F.relu(self.fc46(x)))
        o7 = self.fc57(F.relu(self.fc47(x)))
        o8 = self.fc58(F.relu(self.fc48(x)))
        o9 = self.fc59(F.relu(self.fc49(x)))
        o10 = self.fc510(F.relu(self.fc410(x)))
        o11 = self.fc511(F.relu(self.fc411(x)))
        return o0,o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,o11

"""
Structure of the 24 MB-DNN for cTPnet.
"""
class Net24(nn.Module):
    def __init__(self):
        super(Net24, self).__init__()
        self.fc1 = nn.Linear(12363, 1000)
        self.fc3 = nn.Linear(1000, 128)
        self.fc40 = nn.Linear(128, 64)
        self.fc50 = nn.Linear(64, 1)
        self.fc41 = nn.Linear(128, 64)
        self.fc51 = nn.Linear(64, 1)
        self.fc42 = nn.Linear(128, 64)
        self.fc52 = nn.Linear(64, 1)
        self.fc43 = nn.Linear(128, 64)
        self.fc53 = nn.Linear(64, 1)
        self.fc44 = nn.Linear(128, 64)
        self.fc54 = nn.Linear(64, 1)
        self.fc45 = nn.Linear(128, 64)
        self.fc55 = nn.Linear(64, 1)
        self.fc46 = nn.Linear(128, 64)
        self.fc56 = nn.Linear(64, 1)
        self.fc47 = nn.Linear(128, 64)
        self.fc57 = nn.Linear(64, 1)
        self.fc48 = nn.Linear(128, 64)
        self.fc58 = nn.Linear(64, 1)
        self.fc49 = nn.Linear(128, 64)
        self.fc59 = nn.Linear(64, 1)
        self.fc410 = nn.Linear(128, 64)
        self.fc510 = nn.Linear(64, 1)
        self.fc411 = nn.Linear(128, 64)
        self.fc511 = nn.Linear(64, 1)
        self.fc412 = nn.Linear(128, 64)
        self.fc512 = nn.Linear(64, 1)
        self.fc413 = nn.Linear(128, 64)
        self.fc513 = nn.Linear(64, 1)
        self.fc414 = nn.Linear(128, 64)
        self.fc514 = nn.Linear(64, 1)
        self.fc415 = nn.Linear(128, 64)
        self.fc515 = nn.Linear(64, 1)
        self.fc416 = nn.Linear(128, 64)
        self.fc516 = nn.Linear(64, 1)
        self.fc417 = nn.Linear(128, 64)
        self.fc517 = nn.Linear(64, 1)
        self.fc418 = nn.Linear(128, 64)
        self.fc518 = nn.Linear(64, 1)
        self.fc419 = nn.Linear(128, 64)
        self.fc519 = nn.Linear(64, 1)
        self.fc420 = nn.Linear(128, 64)
        self.fc520 = nn.Linear(64, 1)
        self.fc421 = nn.Linear(128, 64)
        self.fc521 = nn.Linear(64, 1)
        self.fc422 = nn.Linear(128, 64)
        self.fc522 = nn.Linear(64, 1)
        self.fc423 = nn.Linear(128, 64)
        self.fc523 = nn.Linear(64, 1)
    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc3(x))
        o0 = self.fc50(F.relu(self.fc40(x)))
        o1 = self.fc51(F.relu(self.fc41(x)))
        o2 = self.fc52(F.relu(self.fc42(x)))
        o3 = self.fc53(F.relu(self.fc43(x)))
        o4 = self.fc54(F.relu(self.fc44(x)))
        o5 = self.fc55(F.relu(self.fc45(x)))
        o6 = self.fc56(F.relu(self.fc46(x)))
        o7 = self.fc57(F.relu(self.fc47(x)))
        o8 = self.fc58(F.relu(self.fc48(x)))
        o9 = self.fc59(F.relu(self.fc49(x)))
        o10 = self.fc510(F.relu(self.fc410(x)))
        o11 = self.fc511(F.relu(self.fc411(x)))
        o12 = self.fc512(F.relu(self.fc412(x)))
        o13 = self.fc513(F.relu(self.fc413(x)))
        o14 = self.fc514(F.relu(self.fc414(x)))
        o15 = self.fc515(F.relu(self.fc415(x)))
        o16 = self.fc516(F.relu(self.fc416(x)))
        o17 = self.fc517(F.relu(self.fc417(x)))
        o18 = self.fc518(F.relu(self.fc418(x)))
        o19 = self.fc519(F.relu(self.fc419(x)))
        o20 = self.fc520(F.relu(self.fc420(x)))
        o21 = self.fc521(F.relu(self.fc421(x)))
        o22 = self.fc522(F.relu(self.fc422(x)))
        o23 = self.fc523(F.relu(self.fc423(x)))
        return o0,o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,o11,o12,o13,o14,o15,o16,o17,o18,o19,o20,o21,o22,o23




"""
Structure of the n protein class for custom prediction.
"""

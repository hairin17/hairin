subs = [7
9
11
12
13
17
18
19
20
21
27
28
29
30
34
35
36
37
38
39
41
42
45
46
50
53
54
55
56
57
60
61
62
63
64
65
69
70
73
76
77
78
79
80
82
83
84
87
88
89
90
91
92
95
96
97
98
101
104
105
106
108
109
110
112
113
114
115
116
117
119
120
123
124
125]

    for s = 1:numel(subs);
    subj = subs(s)
    dirname = [num2str(subj,'%03d')]
    mkdir(dirname);
    end
    
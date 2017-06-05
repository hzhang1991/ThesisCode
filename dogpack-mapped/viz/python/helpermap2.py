def average(px1,px2,px3,px4,py1,py2,py3,py4,a):
    ave=0.0;
    ave =  ((-2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * py2 * py2 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * py3 * py3 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * py3 * py3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * py3 * py3 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * py3 * py3 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * py3 * py3 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * py3 * py3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * py4 * py4 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * py4 * py4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * py4 * py4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * py4 * py4 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * py4 * py4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * py4 * py4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px1 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px1 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px1 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px1 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px2 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px2 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px2 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px4 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px4 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * py1 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * py1 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * py1 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * py1 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * py1 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * py1 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * py2 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * py2 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * py2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * py2 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * py2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * py2 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * py4 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * py4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * py4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * py4 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * py4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * py4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px2 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px2 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px2 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px2 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px3 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px3 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px3 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px4 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px4 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * py2 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * py2 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * py2 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * py2 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * py2 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * py2 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * py3 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * py3 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * py3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * py3 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * py3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * py3 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * py4 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * py4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * py4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * py4 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * py4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * py4 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px1 * px1 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * px1 * px1 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px1 * px1 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px1 * px1 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * px1 * px1 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px1 * px1 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px2 * px2 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * px2 * px2 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px2 * px2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px2 * px2 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * px2 * px2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px2 * px2 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px4 * px4 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * px4 * px4 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px4 * px4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px4 * px4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * px4 * px4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px4 * px4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * py1 * py1 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * py1 * py1 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * py1 * py1 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * py1 * py1 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * py1 * py1 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * py1 * py1 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * py2 * py2 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * py2 * py2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * py2 * py2 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * py2 * py2 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * py2 * py2 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * py2 * py2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * py4 * py4 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * py4 * py4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * py4 * py4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * py4 * py4 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * py4 * py4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * py4 * py4 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px2 * px2 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * px2 * px2 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px2 * px2 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px2 * px2 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * px2 * px2 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px2 * px2 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px3 * px3 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * px3 * px3 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px3 * px3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px3 * px3 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * px3 * px3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px3 * px3 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px4 * px4 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * px4 * px4 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px4 * px4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px4 * px4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * px4 * px4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px4 * px4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * py2 * py2 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * py2 * py2 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * py2 * py2 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * py2 * py2 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * py2 * py2 + 36 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] - 6 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] + 6 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] + 6 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] - 6 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] + 36 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] - 6 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] + 6 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] + 6 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] - 6 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px1 * px2 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * px1 * px2 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px1 * px2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px1 * px2 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * px1 * px2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px1 * px2 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px1 * px4 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * px1 * px4 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px1 * px4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px1 * px4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * px1 * px4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px1 * px4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px1 * py1 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px1 * py1 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px1 * py1 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px1 * py1 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px1 * py2 - abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px1 * py2 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px1 * py2 - abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px1 * py2 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px1 * py4 - abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px1 * py4 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px1 * py4 - abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px1 * py4 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px2 * px4 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * px2 * px4 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px2 * px4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px2 * px4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * px2 * px4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px2 * px4 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px2 * py1 - abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px2 * py1 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px2 * py1 - abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px2 * py1 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px2 * py2 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px2 * py2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px2 * py2 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px2 * py2 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px2 * py4 - abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px2 * py4 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px2 * py4 - abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px2 * py4 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px4 * py1 - abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px4 * py1 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px4 * py1 - abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px4 * py1 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px4 * py2 - abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px4 * py2 + abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px4 * py2 - abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px4 * py2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * px4 * py4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * px4 * py4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * px4 * py4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * px4 * py4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * py1 * py2 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * py1 * py2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * py1 * py2 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * py1 * py2 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * py1 * py2 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * py1 * py2 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * py1 * py4 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * py1 * py4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * py1 * py4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * py1 * py4 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * py1 * py4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * py1 * py4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[0] * py2 * py4 - 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[1] * py2 * py4 + 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[2] * py2 * py4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[3] * py2 * py4 + 4 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[4] * py2 * py4 - 2 * abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) * a[5] * py2 * py4 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px2 * px3 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * px2 * px3 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px2 * px3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px2 * px3 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * px2 * px3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px2 * px3 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px2 * px4 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * px2 * px4 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px2 * px4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px2 * px4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * px2 * px4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px2 * px4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px2 * py2 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px2 * py2 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px2 * py2 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px2 * py2 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px2 * py3 - abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px2 * py3 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px2 * py3 - abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px2 * py3 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px2 * py4 - abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px2 * py4 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px2 * py4 - abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px2 * py4 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px3 * px4 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * px3 * px4 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px3 * px4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px3 * px4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * px3 * px4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px3 * px4 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px3 * py2 - abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px3 * py2 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px3 * py2 - abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px3 * py2 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px3 * py3 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px3 * py3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px3 * py3 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px3 * py3 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px3 * py4 - abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px3 * py4 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px3 * py4 - abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px3 * py4 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px4 * py2 - abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px4 * py2 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px4 * py2 - abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px4 * py2 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px4 * py3 - abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px4 * py3 + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px4 * py3 - abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px4 * py3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * px4 * py4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * px4 * py4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * px4 * py4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * px4 * py4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * py2 * py3 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * py2 * py3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * py2 * py3 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * py2 * py3 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * py2 * py3 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * py2 * py3 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * py2 * py4 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * py2 * py4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * py2 * py4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * py2 * py4 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * py2 * py4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * py2 * py4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[0] * py3 * py4 - 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[1] * py3 * py4 + 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[2] * py3 * py4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[3] * py3 * py4 + 4 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[4] * py3 * py4 - 2 * abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) * a[5] * py3 * py4) / (abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) + abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3))) / 0.36e2;
    almostarea= abs(px1 * py2 - px1 * py4 - px2 * py1 + px2 * py4 + px4 * py1 - px4 * py2) / 0.2e1 +  abs(px2 * py3 - px2 * py4 - px3 * py2 + px3 * py4 + px4 * py2 - px4 * py3) / 0.2e1;
    return ave;

#----------------------------------------------------------
def get_grid_type(outputdir):

    import string

    Fparams = "".join((outputdir,"/qhelp.dat"     ))
    Rparams = open(Fparams,'r')

    linestring = Rparams.readline()
    linelist = string.split(linestring)
    ndims = int(linelist[0])

    if ndims!=2:
        print ""
        print " Incorrect dimension, ndims must be 2. ndims = ",ndims
        print ""
        return -1

    linestring = Rparams.readline()
    linelist = string.split(linestring)
    GridType = str(linelist[0])

    Rparams.close()

    return GridType
#----------------------------------------------------------
    

#----------------------------------------------------------
def read_params(outputdir,params):

    import string

    Fparams = "".join((outputdir,"/qhelp.dat"     ))
    Rparams = open(Fparams,'r')

    linestring = Rparams.readline()
    linelist = string.split(linestring)
    ndims = int(linelist[0])

    linestring = Rparams.readline()
    linelist = string.split(linestring)
    GridType = str(linelist[0])

    if (GridType=="Cartesian"):
        for k in range (0,11):
            linestring = Rparams.readline()
            linelist = string.split(linestring)
            params[k] = float(linelist[0])
    elif (GridType=="Unstructured"):
        for k in range (0,12):
            linestring = Rparams.readline()
            linelist = string.split(linestring)
            params[k] = int(linelist[0])

    Rparams.close()

    return ndims
#----------------------------------------------------------


#----------------------------------------------------------
def get_kmax(meth1,ndims):

    if (ndims==1):
        return meth1
    elif (ndims==2):
        return int((meth1*(meth1+1))/2)
    elif (ndims==3):
        return int((meth1*(meth1+1)*(meth1+2))/6)
    else:
        print ""
        print " Incorrect dimension in get_kmax, ndims must be 1, 2, or 3. ndims = ",ndims
        print ""
        return -1
    
#----------------------------------------------------------


#----------------------------------------------------------
def read_qfile(mtmp,qfile,qtmp):

    import string

    # open file
    Rqfile = open(qfile,'r')

    # get time
    linestring = Rqfile.readline()
    linelist = string.split(linestring)
    time = float(linelist[0])
    
    # store all Legendre coefficients in qtmp
    for k in range (0,mtmp):
        linestring = Rqfile.readline()
        linelist = string.split(linestring)
        qtmp[k] = float(linelist[0])

    # close file
    Rqfile.close()

    # return time
    return time
#----------------------------------------------------------


#----------------------------------------------------------
#  Sample Legendre polynomial on the midpoint of each element
def GetCart2Legendre(meth1,points_per_dir,s2d,LegVals):

    from math import sqrt
  
    sq3 = sqrt(3.0)
    sq5 = sqrt(5.0)
    sq7 = sqrt(7.0)

    for m in range(0,points_per_dir*points_per_dir):
      xi  = s2d[m,0]
      eta = s2d[m,1]

      xi2 = xi*xi
      xi3 = xi2*xi
      xi4 = xi3*xi

      eta2 = eta*eta
      eta3 = eta2*eta
      eta4 = eta3*eta

      if (meth1==5):
        LegVals[0,m]  = 1.0      
        LegVals[1,m]  = sq3*xi
        LegVals[2,m]  = sq3*eta      
        LegVals[3,m]  = 3.0*xi*eta 
        LegVals[4,m]  = sq5*(1.5*xi2 - 0.5)
        LegVals[5,m]  = sq5*(1.5*eta2 - 0.5)      
        LegVals[6,m]  = sq3*sq5*eta*(1.5*xi2 - 0.5)      
        LegVals[7,m]  = sq3*sq5*xi*(1.5*eta2 - 0.5)      
        LegVals[8,m]  = sq7*(2.5*xi3 - 1.5*xi)
        LegVals[9,m]  = sq7*(2.5*eta3 - 1.5*eta)
        LegVals[10,m] = sq3*sq7*(2.5*xi3 - 1.5*xi)*eta
        LegVals[11,m] = sq3*sq7*(2.5*eta3 - 1.5*eta)*xi
        LegVals[12,m] = 5.0/4.0*(3.0*xi2 - 1.0)*(3.0*eta2 - 1.0)      
        LegVals[13,m] = 105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0      
        LegVals[14,m] = 105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0

      elif (meth1==4):
        LegVals[0,m]  = 1.0      
        LegVals[1,m]  = sq3*xi
        LegVals[2,m]  = sq3*eta      
        LegVals[3,m]  = 3.0*xi*eta 
        LegVals[4,m]  = sq5*(1.5*xi2 - 0.5)
        LegVals[5,m]  = sq5*(1.5*eta2 - 0.5)      
        LegVals[6,m]  = sq3*sq5*eta*(1.5*xi2 - 0.5)      
        LegVals[7,m]  = sq3*sq5*xi*(1.5*eta2 - 0.5)      
        LegVals[8,m]  = sq7*(2.5*xi3 - 1.5*xi)
        LegVals[9,m]  = sq7*(2.5*eta3 - 1.5*eta)

      elif (meth1==3):
        LegVals[0,m]  = 1.0
        LegVals[1,m]  = sq3*xi
        LegVals[2,m]  = sq3*eta
        LegVals[3,m]  = 3.0*xi*eta
        LegVals[4,m]  = sq5*(1.5*xi2 - 0.5)
        LegVals[5,m]  = sq5*(1.5*eta2 - 0.5)

      elif (meth1==2):
        LegVals[0,m]  = 1.0  
        LegVals[1,m]  = sq3*xi
        LegVals[2,m]  = sq3*eta

      elif (meth1==1):
        LegVals[0,m]  = 1.0
#----------------------------------------------------------


#----------------------------------------------------------
# Turn coefficients into point values on 2D Cartesian grid
#
"""
def sample_state2_cart_mod(mx_old,my_old,points_per_dir,meqn,kmax,qcoeffs,LegVals,qsoln):

    index = 0

    for j in range(1,my_old+1):
        for m1 in range(1,points_per_dir+1):
            for i in range(1,mx_old+1):
                for m2 in range(1,points_per_dir+1):
    
                    m = m2 + points_per_dir*(m1-1)
                    index = index + 1
    
                    for n in range(1,meqn+1):
                        qsoln[index-1,n-1] = 0.0

                        for k in range(1,kmax+1):
                            qsoln[index-1,n-1] = qsoln[index-1,n-1] + qcoeffs[k-1,n-1,j-1,i-1]*LegVals[k-1,m-1]
   
#----------------------------------------------------------
"""

#----------------------------------------------------------
# Turn coefficients into point values on 2D Cartesian grid
#
def sample_state2_cart_mod_map(mx_old,my_old,points_per_dir,meqn,kmax,qcoeffs,LegVals,qsoln,xpl,ypl):
    import numpy
    index = 0
    a=numpy.zeros((6,));
    for j in range(1,my_old+1):
        for m1 in range(1,points_per_dir+1):
            for i in range(1,mx_old+1):
                for m2 in range(1,points_per_dir+1):
    
                    m = m2 + points_per_dir*(m1-1)
                    index = index + 1
    
                    for n in range(1,meqn+1):
                       
                        xp1 = xpl[i-1,j-1];
                        yp1=  ypl[i-1,j-1];
                        xp2 = xpl[i,j-1];
                        yp2=  ypl[i,j-1];
                        xp3 = xpl[i,j];
                        yp3=  ypl[i,j];
                        xp4 = xpl[i-1,j];
                        yp4=  ypl[i-1,j];
                        qsoln[index-1,n-1] = 0.0
                        a[0]=qcoeffs[0,n-1,j-1,i-1];
                        a[1]=qcoeffs[1,n-1,j-1,i-1];
                        a[2]=qcoeffs[2,n-1,j-1,i-1];
                        a[3]=qcoeffs[3,n-1,j-1,i-1];
                        a[4]=qcoeffs[4,n-1,j-1,i-1];
                        a[5]=qcoeffs[5,n-1,j-1,i-1];
                        area=(0.1250000000*yp4*xp3-0.1250000000*xp4*yp3-0.1250000000*yp4*xp1+0.1250000000*xp4*yp1+0.1250000000*yp3*xp2-0.1250000000*yp1*xp2-0.1250000000*xp3*yp2+0.1250000000*xp1*yp2)*4.0
                        qsoln[index-1,n-1] = area*a[0]#average(xp1,xp2,xp3,xp4,yp1,yp2,yp3,yp4,a)

def sample_state2_cart_mod(mx_old,my_old,points_per_dir,meqn,kmax,qcoeffs,LegVals,qsoln):

    index = 0

    for j in range(1,my_old+1):
        for m1 in range(1,points_per_dir+1):
            for i in range(1,mx_old+1):
                for m2 in range(1,points_per_dir+1):
    
                    m = m2 + points_per_dir*(m1-1)
                    index = index + 1
    
                    for n in range(1,meqn+1):
                        qsoln[index-1,n-1] = 0.0

                        for k in range(1,kmax+1):
                            qsoln[index-1,n-1] = qsoln[index-1,n-1] + qcoeffs[k-1,n-1,j-1,i-1]*LegVals[k-1,m-1]

#----------------------------------------------------------
def read_mesh_params(meshdir,mesh_params):

    import string

    Fparams = "".join((meshdir,"/mesh_params.dat"))
    Rparams = open(Fparams,'r')

    for k in range (0,7):
        linestring = Rparams.readline()
        linelist = string.split(linestring)
        mesh_params[k] = int(linelist[0])

    Rparams.close()

#----------------------------------------------------------


#----------------------------------------------------------
def read_node(meshdir,NumPhysNodes,x,y):  

    import string

    Fnode   = "".join((meshdir,"/mesh_node.dat"       ))
    Rnode   = open(Fnode,  'r')

    for i in range(0,NumPhysNodes):
        linestring = Rnode.readline()
        linelist = string.split(linestring)
        x[i] = float(linelist[0])
        y[i] = float(linelist[1])
        
    Rnode.close()
#----------------------------------------------------------


#----------------------------------------------------------
def read_tnode(meshdir,NumPhysElems,tnode):

    import string
    
    Ftnode  = "".join((meshdir,"/mesh_tnode.dat"      ))
    Rtnode  = open(Ftnode, 'r')
    
    for i in range(0,NumPhysElems):
        linestring = Rtnode.readline()
        linelist = string.split(linestring)
        tnode[i,0] = int(linelist[0])-1
        tnode[i,1] = int(linelist[1])-1
        tnode[i,2] = int(linelist[2])-1

    Rtnode.close()
#----------------------------------------------------------


#----------------------------------------------------------
def set_tcounter(NumPhysElems,tnode,tcounter):

    for i in range(0,NumPhysElems):
        j1 = tnode[i,0]
        j2 = tnode[i,1]
        j3 = tnode[i,2]

        tcounter[j1] = tcounter[j1] + 1
        tcounter[j2] = tcounter[j2] + 1
        tcounter[j3] = tcounter[j3] + 1

#----------------------------------------------------------
    

#----------------------------------------------------------
def GetMonomialToLegendre(kmax,Mon2Leg):

    from math import sqrt
    
    sq2 = sqrt(2.0);
    sq3 = sqrt(3.0);
    sq5 = sqrt(5.0);
    sq7 = sqrt(7.0);
    sq10 = sqrt(10.0);
    sq13 = sqrt(13.0);
    sq19 = sqrt(19.0);
    sq23 = sqrt(23.0);
    sq71 = sqrt(71.0);

    Mmat = [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            #
            [0.0, 3.0*sq2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            #
            [0.0, sq2*sq3, 2.0*sq2*sq3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            #
            [(5.0/21.0)*sq7, 8.0/sq7, 8.0/sq7, 60.0/sq7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            #
            [-10.0/9.0*(sq3/sq7), -2.0/(sq3*sq7), 4.0*sq3/sq7, 30.0*sq3/sq7, 5.0*sq3*sq7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            #
            [-2.0/9.0*sq3*sq5, 2.0*sq5/sq3, 0, 6.0*sq3*sq5, sq3*sq5, 6.0*sq3*sq5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            #
            [-8.0*sq3*sq13/117.0, 20.0/(sq3*sq13), -20.0/(sq3*sq13), 20.0*sq3/sq13, 40.0*sq3/sq13, 0.0, 210.0*sq3/sq13, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            #
            [-8.0/351.0*sq3*sq5*sq13*sq19, -4.0/247.0*sq3*sq5*sq13*sq19, 4.0/247.0*sq3*sq5*sq13*sq19, (20.0/39.0)*sq3*sq5*sq13*sq19, (344.0/741.0)*sq3*sq5*sq13*sq19, (32.0/57.0)*sq3*sq5*sq13*sq19, (602.0/247.0)*sq3*sq5*sq13*sq19, (56.0/19.0)*sq3*sq5*sq13*sq19, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            #
            [0.0, -52.0/(sq3*sq19), -(8.0/19.0)*sq3*sq19, (8.0/3.0)*sq3*sq19, -4.0/(sq3*sq19), 80.0/(sq3*sq19), (392.0/19.0)*sq3*sq19, (140.0/19.0)*sq3*sq19, 14.0*sq3*sq19, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            #
            [0.0, -(4.0/3.0)*sq7, -(8.0/3.0)*sq7, 8.0*sq7, 4.0*sq7, 0, 24.0*sq7, 60.0*sq7, 2.0*sq7, 40.0*sq7, 0.0, 0.0, 0.0, 0.0, 0.0],
            #
            [-(2.0/27.0)*sq5*sq7, -(26.0/21.0)*sq5*sq7, -(4.0/21.0)*sq5*sq7, -6.0*sq5*sq7, sq5*sq7, 0.0, 0.0, 0.0, 16.0*sq5*sq7, 0, 72.0*sq5*sq7, 0.0, 0.0, 0.0, 0.0],
            #
            [-(2.0/1053.0)*sq5*sq7*sq13*sq23, (314.0/6279.0)*sq5*sq7*sq13*sq23, -(544.0/6279.0)*sq5*sq7*sq13*sq23, -(2.0/13.0)*sq5*sq7*sq13*sq23, -(47.0/897.0)*sq5*sq7*sq13*sq23, (70.0/897.0)*sq5*sq7*sq13*sq23, 0.0, 0.0, -(752.0/897.0)*sq5*sq7*sq13*sq23, (1120.0/897.0)*sq5*sq7*sq13*sq23, -(1128.0/299.0)*sq5*sq7*sq13*sq23, (1680.0/299.0)*sq5*sq7*sq13*sq23, 0.0, 0.0, 0.0],
            #
            [(55.0/74763.0)*sq7*sq13*sq71, -(440.0/2769.0)*sq7*sq13*sq71, -(440.0/2769.0)*sq7*sq13*sq71, -(140.0/923.0)*sq7*sq13*sq71, -(160.0/2769.0)*sq7*sq13*sq71, -(160.0/2769.0)*sq7*sq13*sq71, (360.0/71.0)*sq7*sq13*sq71, (360.0/71.0)*sq7*sq13*sq71, (6800.0/2769.0)*sq7*sq13*sq71, (6800.0/2769.0)*sq7*sq13*sq71, (10200.0/923.0)*sq7*sq13*sq71, (10200.0/923.0)*sq7*sq13*sq71, (1620.0/71.0)*sq7*sq13*sq71, 0.0, 0.0],
            #
            [(2.0/639.0)*sq3*sq5*sq23*sq71, (4.0/14697.0)*sq3*sq5*sq23*sq71, -(280.0/14697.0)*sq3*sq5*sq23*sq71, -(20.0/71.0)*sq3*sq5*sq23*sq71, -(466.0/1633.0)*sq3*sq5*sq23*sq71, -(40.0/1633.0)*sq3*sq5*sq23*sq71, (60.0/71.0)*sq3*sq5*sq23*sq71, (60.0/sq71)*sq3*sq5*sq23, -(4.0/1633.0)*sq3*sq5*sq23*sq71, (280.0/1633.0)*sq3*sq5*sq23*sq71, (9780.0/1633.0)*sq3*sq5*sq23*sq71, (1260.0/1633.0)*sq3*sq5*sq23*sq71, (270.0/sq71)*sq3*sq5*sq23, 3.0*sq3*sq5*sq23*sq71, 0.0],
            #
            [(2.0/9.0)*sq5, -(4.0/3.0)*sq5, 0.0, -20.0*sq5, -2.0*sq5, -20.0*sq5, 60.0*sq5, 60.0*sq5, 12.0*sq5, 0.0, 60.0*sq5, 420.0*sq5, 270.0*sq5, 3.0*sq5, 210.0*sq5]]

    for i in range(0,kmax):
        for j in range(0,kmax):
            Mon2Leg[i,j] = Mmat[i][j]

#----------------------------------------------------------


#----------------------------------------------------------
def GetUnstLegendre(meth1,kmax,points_per_dir,zx,zy,Mon2Leg,MonVals,LegVals):

    p2 = points_per_dir*points_per_dir
    
    for m in range (0,p2):
        xi  = zx[m]
        eta = zy[m]

        if (meth1>0):
            MonVals[0,m]  = 1.0;

        if (meth1>1):
            MonVals[1,m]  = xi;
            MonVals[2,m]  = eta;

        if (meth1>2):
            MonVals[3,m]  = xi*eta;
            MonVals[4,m]  = xi*xi;
            MonVals[5,m]  = eta*eta;

        if (meth1>3):
            MonVals[6,m]  = eta*xi*xi;
            MonVals[7,m]  = xi*eta*eta;
            MonVals[8,m]  = xi*xi*xi;
            MonVals[9,m]  = eta*eta*eta;

        if (meth1>4):
            MonVals[10,m] = xi*xi*xi*eta;
            MonVals[11,m] = eta*eta*eta*xi;
            MonVals[12,m] = xi*xi*eta*eta;
            MonVals[13,m] = xi*xi*xi*xi;
            MonVals[14,m] = eta*eta*eta*eta;
    
    for m in range(0,p2):
        for k in range(0,kmax):
            LegVals[k,m] = 0.0
            for kk in range(0,kmax):
                LegVals[k,m] = LegVals[k,m] + Mon2Leg[k,kk]*MonVals[kk,m]
#----------------------------------------------------------


#----------------------------------------------------------
def set_soln_at_node_values(NumPhysElems,NumPhysNodes,meqn,tcounter,tnode,qsoln_elem,qsoln):

    for i in range(0,NumPhysElems):
        j1 = tnode[i,0]
        j2 = tnode[i,1]
        j3 = tnode[i,2]

        for m in range(0,meqn):
            qsoln[j1,m] = qsoln[j1,m] + qsoln_elem[i,m]
            qsoln[j2,m] = qsoln[j2,m] + qsoln_elem[i,m]
            qsoln[j3,m] = qsoln[j3,m] + qsoln_elem[i,m]

    for i in range(0,NumPhysNodes):
        for m in range(0,meqn):
            qsoln[i,m] = qsoln[i,m]/float(tcounter[i])

#----------------------------------------------------------


#----------------------------------------------------------
def sample_state2_unst(NumPhysElems,mpts,meqn,kmax,LegVals,qcoeffs,qsoln_elem):

    index = -1

    for i in range(0,NumPhysElems):
        for m in range(0,mpts):
            index = index+1
            for j in range(0,meqn):
                for k in range(0,kmax):
                    qsoln_elem[index,j] = qsoln_elem[index,j] + qcoeffs[i,j,k]*LegVals[k,m]

#----------------------------------------------------------


#----------------------------------------------------------
def DivideUnstMesh(points_per_dir,NumPhysElems,x,y,tnode,
                   x_new,y_new,tnode_new,zx,zy,newsizes):

    # 
    #  Add extra points and elements if points_per_dir>1
    #

    import numpy as np
    from math import pi
    from math import atan2
    from math import sqrt
    from math import fabs

    onethird = 1.0/3.0
    twothird = 2.0/3.0
    points_per_elem = ((points_per_dir+1)*(points_per_dir+2))/2
    node_list = np.zeros(NumPhysElems*points_per_elem,float)
    map_list  = np.zeros(NumPhysElems*points_per_elem,int)

    h = 1.0/float(points_per_dir)
    is_empty_flag = True
    node_list_size = -1   # new node index
    map_counter = -1      # helper array index
    elem_list_size = -1   # new element index

    tmp_stopper = 0
    k = 0
    
    for n in range (0,NumPhysElems):
      
        x1 = x[tnode[n,0]]
        y1 = y[tnode[n,0]]
      
        x2 = x[tnode[n,1]]
        y2 = y[tnode[n,1]]
      
        x3 = x[tnode[n,2]]
        y3 = y[tnode[n,2]]
      
        xm = onethird*(x1+x2+x3)
        ym = onethird*(y1+y2+y3)
      
        km = np.zeros(points_per_elem,int)
        index = -1
      
        for j in range(0,points_per_dir+1):
            for i in range(0,(points_per_dir+1-j)):
                index = index+1;      
                k = k+1
                km[index] = k
          
                xi  = float(i)*h-onethird
                eta = float(j)*h-onethird
          
                xtmp = xm + (x2-x1)*xi + (x3-x1)*eta
                ytmp = ym + (y2-y1)*xi + (y3-y1)*eta
          
                ztmp = 3*pi + atan2(ytmp,xtmp) + 1.0e3*sqrt(pow(xtmp,2)+pow(ytmp,2));
          
                if (is_empty_flag):
                    kk = 0
                    node_list[0] = ztmp
                    node_list_size = node_list_size + 1
                    map_counter = map_counter+1
                    map_list[map_counter]  = kk
            
                    x_new[kk] = xtmp
                    y_new[kk] = ytmp

                    is_empty_flag = False
                else:
                    map_counter = map_counter+1
                    min_value = fabs(ztmp-node_list[0])
                    min_value_location = 0
                    mstop = 0
                    ijk = 0
                    while(mstop==0 and node_list_size>1):
                        ijk = ijk+1
                        tmp_value = fabs(ztmp-node_list[ijk])
                        if tmp_value<min_value:
                            min_value = tmp_value
                            min_value_location = ijk
                        if (tmp_value<=1.0e-8 or ijk==node_list_size):
                            mstop = 1

                    if ( min_value >= 1.0e-8 ):
                        kk = kk+1
                        node_list[kk] = ztmp
                        node_list_size = node_list_size + 1
                        
                        map_list[map_counter] = kk
                        x_new[kk] = xtmp
                        y_new[kk] = ytmp
                    else:
                        map_list[map_counter] = min_value_location
                        
        # all triangles oriented one way
        m1 = km[0]
        m2 = km[points_per_dir+1]
        for j in range(1,points_per_dir+1):
            for i in range(1,(points_per_dir+2-j)):

                elem_list_size = elem_list_size + 1
                tnode_new[elem_list_size,0] = map_list[m1+i-2]
                tnode_new[elem_list_size,1] = map_list[m1+i-1]
                tnode_new[elem_list_size,2] = map_list[m2+i-2]
        
            m1 = m1 + (points_per_dir+2-j)
            m2 = m2 + (points_per_dir+1-j)

        # all triangles oriented the other way
        m1 = km[0]
        m2 = km[points_per_dir+1]
        for j in range(1,points_per_dir):
            for i in range(1,(points_per_dir+1-j)):

                elem_list_size = elem_list_size + 1
                tnode_new[elem_list_size,0] = map_list[m1+i-1]
                tnode_new[elem_list_size,1] = map_list[m2+i-1]
                tnode_new[elem_list_size,2] = map_list[m2+i-2]
        
            m1 = m1 + (points_per_dir+2-j)
            m2 = m2 + (points_per_dir+1-j)

    NumPhysElems_new = elem_list_size+1
    NumPhysNodes_new = node_list_size+1

    newsizes[0] = NumPhysElems_new
    newsizes[1] = NumPhysNodes_new

    # Properly orient each triangle
    for j in range(0,NumPhysElems_new+1):

        x1 = x_new[tnode_new[n,0]]
        y1 = y_new[tnode_new[n,0]]
      
        x2 = x_new[tnode_new[n,1]]
        y2 = y_new[tnode_new[n,1]]
      
        x3 = x_new[tnode_new[n,2]]
        y3 = y_new[tnode_new[n,2]]

        Area = (x3-x2)*y1 + (x1-x3)*y2 + (x2-x1)*y3

        if (Area<0.0):
            tmp = tnode_new[j,1]
            tnode_new[j,1] = tnode_new[j,2]
            tnode_new[j,2] = tmp

    # Set sample point values on each element
    x1 = x[tnode[0,0]]
    y1 = y[tnode[0,0]]
    
    x2 = x[tnode[0,1]]
    y2 = y[tnode[0,1]]
    
    x3 = x[tnode[0,2]]
    y3 = y[tnode[0,2]]
    
    Area = (x3-x2)*y1 + (x1-x3)*y2 + (x2-x1)*y3

    for n in range(0,points_per_dir*points_per_dir):

        p1 = onethird*(x_new[tnode_new[n,0]]+x_new[tnode_new[n,1]]+x_new[tnode_new[n,2]])
        p2 = onethird*(y_new[tnode_new[n,0]]+y_new[tnode_new[n,1]]+y_new[tnode_new[n,2]])
      
        zx[n] = ((x1-x3)*p2+(-y1+y3)*p1+twothird*y1*x3+onethird*y1*x2-
                 onethird*y2*x1-onethird*y3*x2-twothird*y3*x1+onethird*y2*x3)/Area
      
        zy[n] = ((-x1+x2)*p2+(y1-y2)*p1+twothird*y2*x1+onethird*y3*x1-
                 twothird*y1*x2-onethird*y1*x3+onethird*y2*x3-onethird*y3*x2)/Area
 
#----------------------------------------------------------

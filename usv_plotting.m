% plotting resulting data for AI sorting

T = readtable('usv_ex.xls');


duration = T{:,4};
maxfreq = T{:,5};
name = T{1:200,6};

f=1;
fd=1;
c=1;
c2=1;
c3=1;
c4=1;
c5=1;
d=1;
u=1;
h=1;
s=1;
flat = 0;
flat_down = 0;
complex = 0;
complex_2 = 0;
complex_3 = 0;
complex_4 = 0;
complex_5 = 0;
down = 0;
unclassfied = 0;
harmonic = 0;
short = 0;

for i = 1:size(name)

    if name(i) == "flat"
        flat = flat+1;
    end
    if name(i) == "flat down"
        flat_down = flat_down+1;
    end
    if name(i) == "complex"
        complex = complex+1;
    end
    if name(i) == "complex 2"
        complex_2 = complex_2+1;
    end
    if name(i) == "complex 3"
        complex_3 = complex_3+1;
    end
    if name(i) == "complex 4"
        complex_4 = complex_4+1;
    end
    if name(i) == "complex 5"
        complex_5 = complex_5+1;
    end
    if name(i) == "down"
        down = down+1;
    end
    if name(i) == "unclassfied"
        unclassfied = unclassfied+1;
    end
    if name(i) == "harmonic"
        harmonic = harmonic+1;
    end
    if name(i) == "short"
        short = short+1;
    end

end

flat_array = zeros(2,flat);
flat_down_array = zeros(2,flat_down);
complex_array = zeros(2,complex);
complex2_array = zeros(2,complex_2);
complex3_array = zeros(2,complex_3);
complex4_array = zeros(2,complex_4);
complex5_array = zeros(2,complex_5);
down_array = zeros(2,down);
harmonic_array = zeros(2,harmonic);
unclassfied_array = zeros(2,unclassfied);
short_array = zeros(2, short);

for i = 1:size(name)

    if name(i) == "flat"
        flat_array(f,1) = duration(i);
        flat_array(f,2) = maxfreq(i);
        f = f+1;
    end
    if name(i) == "flat down"
        flat_down_array(fd,1) = duration(i);
        flat_down_array(fd,2) = maxfreq(i);
        fd = fd+1;
    end
    if name(i) == "complex"
        complex_array(c,1) = duration(i);
        complex_array(c,2) = maxfreq(i);
        c = c+1;
    end
    if name(i) == "complex 2"
        complex2_array(c2,1) = duration(i);
        complex2_array(c2,2) = maxfreq(i);
        c2 = c2+1;
    end
    if name(i) == "complex 3"
        complex3_array(c3,1) = duration(i);
        complex3_array(c3,2) = maxfreq(i);
        c3 = c3+1;
    end
    if name(i) == "complex 4"
        complex4_array(c4,1) = duration(i);
        complex4_array(c4,2) = maxfreq(i);
        c4 = c4+1;
    end
    if name(i) == "complex 5"
        complex5_array(c5,1) = duration(i);
        complex5_array(c5,2) = maxfreq(i);
        c5 = c5+1;
    end
    if name(i) == "down"
        down_array(d,1) = duration(i);
        down_array(d,2) = maxfreq(i);
        d = d+1;
    end
    if name(i) == "unclassfied"
        unclassfied_array(u,1) = duration(i);
        unclassfied_array(u,2) = maxfreq(i);
        u = u+1;
    end
    if name(i) == "harmonic"
        harmonic_array(h,1) = duration(i);
        harmonic_array(h,2) = maxfreq(i);
        h = h+1;
    end
    if name(i) == "short"
        short_array(s,1) = duration(i);
        short_array(s,2) = maxfreq(i);
        s = s+1;
    end

end



figure
plot(flat_array(:,1),flat_array(:,2),'r*',complex_array(:,1),complex_array(:,2),'go',complex2_array(:,1),complex2_array(:,2),'g+',complex3_array(:,1),complex3_array(:,2),'g*',flat_down_array(:,1),flat_down_array(:,2),'b*',short_array(:,1),short_array(:,2),'*')
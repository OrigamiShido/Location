function data_dec = signal_2(data,M)
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
data_dec=zeros(1,length(data));
for n=1:length(data)

    if data(n)>32767   %8λ��������з�������ȡֵ��Χ[-128,127],(2^8/2-1=127),����ʵ��λ���޸�
       data_dec(n)=double(data(n))-65536;
    else
        data_dec(n) = data(n);
    end

end

end


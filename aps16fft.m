function [res]=apowerspectrum(sequence)
len=length(sequence);
na=0;
nc=0;
ng=0;
nt=0;
for i=1:len
if strcmp(sequence(i),'A')||strcmp(sequence(i),'a')
        ua(i)=1;
        uc(i)=0;
        ug(i)=0;
        ut(i)=0;
        na=na+1;
end
    if strcmp(sequence(i),'C')||strcmp(sequence(i),'c')
        ua(i)=0;
        uc(i)=1;
        ug(i)=0;
        ut(i)=0;
        nc=nc+1;
    end
    if strcmp(sequence(i),'G')||strcmp(sequence(i),'g')
        ua(i)=0;
        uc(i)=0;
        ug(i)=1;
        ut(i)=0;
        ng=ng+1;
    end
     if strcmp(sequence(i),'T')||strcmp(sequence(i),'t')
        ua(i)=0;
        uc(i)=0;
        ug(i)=0;
        ut(i)=1;
        nt=nt+1;
     end
end

UA=fft(ua);
UC=fft(uc);
UG=fft(ug);
UT=fft(ut);

PSA=zeros(1,len);
PSC=zeros(1,len);
PSG=zeros(1,len);
PST=zeros(1,len);
for i=1:len
    PSA(1,i)=(abs(UA(i)))^2;
    PSC(1,i)=(abs(UC(i)))^2;
    PSG(1,i)=(abs(UG(i)))^2;
    PST(1,i)=(abs(UT(i)))^2;
end
APSA=zeros(1,len-1);
APSC=zeros(1,len-1);
APSG=zeros(1,len-1);
APST=zeros(1,len-1);
for i=1:(len-1)
    for k=1:i
        APSA(i)=APSA(i)+PSA(k+1);
        APSC(i)=APSC(i)+PSC(k+1);
        APSG(i)=APSG(i)+PSG(k+1);
        APST(i)=APST(i)+PST(k+1);
    end
end

   MA1=0;
   MA2=0;
   AMA1=0;
   AMA2=0;
   MC1=0;
   MC2=0;
   AMC1=0;
   AMC2=0;
   MG1=0;
   MG2=0;
   AMG1=0;
   AMG2=0;
   MT1=0;
   MT2=0;
   AMT1=0;
   AMT2=0;
   for i=1:(len-1)
       MA1=MA1+APSA(i);
       MC1=MC1+APSC(i);
       MG1=MG1+APSG(i);
       MT1=MT1+APST(i);
       MA2=MA2+APSA(i)^2;
       MC2=MC2+APSC(i)^2;
       MG2=MG2+APSG(i)^2;
       MT2=MT2+APST(i)^2;
       AMA1=AMA1+abs(APSA(i)-mean(APSA));
       AMC1=AMC1+abs(APSC(i)-mean(APSC));
       AMG1=AMG1+abs(APSG(i)-mean(APSG));
       AMT1=AMT1+abs(APST(i)-mean(APST));
       AMA2=AMA2+abs(APSA(i)-mean(APSA))^2;
       AMC2=AMC2+abs(APSC(i)-mean(APSC))^2;
       AMG2=AMG2+abs(APSG(i)-mean(APSG))^2;
       AMT2=AMT2+abs(APST(i)-mean(APST))^2;
   end
   MA1=MA1/len;
   MC1=MC1/len;
   MG1=MG1/len;
   MT1=MT1/len;
   
   MA2=MA2/(na*(len-na)*len*len);
   MC2=MC2/(nc*(len-nc)*len*len);
   MG2=MG2/(ng*(len-ng)*len*len);
   MT2=MT2/(nt*(len-nt)*len*len);

   
   AMA1=AMA1/len;
   AMC1=AMC1/len;
   AMG1=AMG1/len;
   AMT1=AMT1/len;
   
   
   AMA2=AMA2/(len*len*na*(len-na));
   AMC2=AMC2/(len*len*nc*(len-nc));
   AMG2=AMG2/(len*len*ng*(len-ng));
   AMT2=AMT2/(len*len*nt*(len-nt));
   
   
   res=[MA1,MA2,AMA1,AMA2,MC1,MC2,AMC1,AMC2,MG1,MG2,AMG1,AMG2,MT1,MT2,AMT1,AMT2];
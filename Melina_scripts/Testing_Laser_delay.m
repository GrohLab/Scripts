% Testing laser delay

% 100x10ms all 2 sec:
%in ADC-00: trigger to activate laser sent by spike 
%in ADC-01: light pulse which is emitted by laser and recorded by transistor 
% see time delay between laser tigger to send out a light pulse and the real light flash coming out of the fiber tip

%%100x10ms all 2 sec with 100% laser intensity
%Datadir='F:\SalienceSetup-Synchronization\Laser testing_delay from trigger
%output and light input\2020-01-14\optoStim_100x10msall2sec\intan'

%%50x10ms al 2 sec 0,10,30,50,70,100
%Datadir='F:\SalienceSetup-Synchronization\Laser testing_delay from trigger output and light input\2020-02-05\0_10_30_50_70_100_50x10msall2sec\intan'

%%50x10ms all2sec 100,70,50,30,10,0
%Datadir='F:\SalienceSetup-Synchronization\Laser testing_delay from trigger output and light input\2020-02-05\100_70_50_30_10_0_50x10msall2sec\intan'


%IF YOU WANT TO LOOK AT ALIGNMENT OF DIG IN AND TRIGGERS
 %figure
 %plot(laserSignal)
 %hold on
 %plot(ltOn,laserSignal(ltOn),'o')

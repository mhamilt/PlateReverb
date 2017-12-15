/*
 ==============================================================================
 
 This file was auto-generated!
 
 ==============================================================================
 */
#include "../../../FDTD_Classes/PlateClass.hpp"
#include "../JuceLibraryCode/JuceHeader.h"

//==============================================================================
/*
 This component lives inside our window, and this is where you should put all
 your controls and content.
 */
class MainContentComponent   :	//public AudioAppComponent,
public Component,
public Slider::Listener,
public AudioIODeviceCallback
{
public:
	
	//==============================================================================
	MainContentComponent() :
	audioSetupComp (audioDeviceManager,
					0,		// minAudioInputChannels
					1,		// maxAudioInputChannels
					0,		// minAudioOutputChannels
					2,		// maxAudioOutputChannels
					false,	// showMidiInputOptions
					false,	// showMidiOutputSelector
					false,	// showChannelsAsStereoPairs
					true),	// hideAdvancedOptionsWithButton
	
	currentSampleRate (0.0)
	
	{
		setSize (600, 300);
		addAndMakeVisible (positionSlider);
		positionSlider.setRange (0.0, 1.0, 0.01);
		positionSlider.setTextBoxStyle (Slider::TextBoxRight, false, 100, 20);
		positionSlider.addListener (this);
		
		addAndMakeVisible (positionLabel);
		positionLabel.setText ("Read Out Point", dontSendNotification);
		
		addAndMakeVisible (levelOutLabel);
		addAndMakeVisible (levelOutSlider);
		levelOutSlider.setRange (0.0, 20.0, 0.01);
		levelOutSlider.setTextBoxStyle (Slider::TextBoxRight, false, 100, 20);
		levelOutSlider.addListener (this);
		levelOutLabel.setText ("Gain Out", dontSendNotification);
		
		audioDeviceManager.initialise (1, 2, 0, true, String(), 0);
		audioDeviceManager.addAudioCallback (this);
		audioDeviceManager.getAudioDeviceSetup(audioIOSettings);
		
		addAndMakeVisible (audioSetupComp);
		currentSampleRate = audioIOSettings.sampleRate;
		plate.setup(currentSampleRate, true);
		plate.setLoss(7,.55);
		currentReadPoint = .1;
		targetReadPoint = currentReadPoint;
		plate.setInterpOut(currentReadPoint,currentReadPoint);
	}
	
	~MainContentComponent()
	{
		//shutdownAudio();
	}
	
	//==============================================================================
	void audioDeviceAboutToStart (AudioIODevice* device) override
	{
		sampleRate = device->getCurrentSampleRate();
		
		inputBuffer.setSize (1, device->getCurrentBufferSizeSamples());
		inputBuffer.clear();
		outputBuffer.setSize (1, device->getCurrentBufferSizeSamples());
		outputBuffer.clear();
	}
	
	
	void audioDeviceIOCallback (const float** inputChannelData, int numInputChannels,
								float** outputChannelData, int numOutputChannels,
								int numSamples) override
	{
		
		float outputSamp;
		if (targetReadPoint != currentReadPoint)
		{
			
			const double smooth_inc = 2e-3	;
			if ((targetReadPoint - currentReadPoint)>smooth_inc)
			{
				currentReadPoint += smooth_inc;
			}
			else if ((currentReadPoint-targetReadPoint)>smooth_inc)
			{
				currentReadPoint -= smooth_inc;
			}
			plate.setInterpOut(currentReadPoint,currentReadPoint);
			
			
			
			for (int i = 0; i < numSamples; ++i)
			{
				
				float inputSamp = 0;
				for (int iSamp = numInputChannels; --iSamp >= 0;)
				{
					if (inputChannelData[iSamp] != 0)
					{
						inputSamp += inputChannelData[iSamp][i];
					}
				}
				outputSamp = plate.reverb(levelOutSliderValue*inputSamp);
				
				for (int j = numOutputChannels; --j >= 0;)
				{
					if (outputChannelData[j] != 0)
					{
						outputChannelData[j][i] = outputSamp;
					}
					
				}
			}
			
		}
		else
		{
			for (int i = 0; i < numSamples; ++i)
			{
				
				float inputSamp = 0;
				float outputSamp;
				for (int iSamp = numInputChannels; --iSamp >= 0;)
				{
					if (inputChannelData[iSamp] != 0)
					{
						inputSamp += inputChannelData[iSamp][i];
					}
				}
				outputSamp = plate.reverb(levelOutSliderValue*inputSamp);
				
				for (int j = numOutputChannels; --j >= 0;)
				{
					if (outputChannelData[j] != 0)
					{
						outputChannelData[j][i] = outputSamp;
					}
					
				}
			}
		}
	}
	
	void audioDeviceStopped() override
	{
	}
	
	
	void sliderValueChanged (Slider* slider) override
	{
		if (slider == &positionSlider)
		{
			targetReadPoint = positionSlider.getValue();
		}
		if (slider == &levelOutSlider)
		{
			levelOutSliderValue = levelOutSlider.getValue();
		}
	}
	
	
	void resized() override
	{
		positionLabel.setBounds (10, 10, 90, 20);
		positionSlider.setBounds (100, 10, getWidth() - 110, 20);
		levelOutLabel.setBounds (10, 30, 90, 20);
		levelOutSlider.setBounds (100, 30, getWidth() - 110, 20);
		audioSetupComp.setBounds(50, 40, 400, 200);
	}
	
private:
	AudioSampleBuffer inputBuffer, outputBuffer;
	AudioDeviceManager audioDeviceManager;
	AudioDeviceSelectorComponent audioSetupComp;
	juce::AudioDeviceManager::AudioDeviceSetup audioIOSettings;
	
	Slider positionSlider, levelOutSlider;
	Label positionLabel, levelOutLabel;
	double currentSampleRate;
	double currentReadPoint, targetReadPoint;
	double levelOutSliderValue = 0.707;
	double readOutIncrement = 0;
	FD_Plate plate;
	double sampleRate;
	
	JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MainContentComponent)
};

//Component* createMainContentComponent()     { return new MainContentComponent(); }


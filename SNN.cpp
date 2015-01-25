// SNN.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdlib.h>
#include <math.h>
#include <time.h>

long refCount = 0;
typedef void(*dispose)(void* args);

void _free(void **pp)
{
	if (pp != NULL)
	{
		free(*pp);
		*pp = NULL;
	}
}

struct eventArgs
{
    void* value;
    dispose dispose;
    long t;
};
typedef void(*delegate)(void*, void*, eventArgs*);

typedef float(*decay)(float potential, long deltaT);
typedef float(*delay)(float signal);
struct neuron
{
    float potential;
    float threshold;
    int refracting;
    long t;
    bool isFiring;
    delegate Integration;
    delegate Activation;
    decay Decay;
    delay Delay;
    struct synapseList *axons; // outputs
    struct synapseList *dendrites; // inputs
};

struct synapse
{
    struct neuron *input;
    struct neuron *output;
    float weight;
};

struct synapseList
{
    synapse *value;
    synapseList *next;
};

struct eventQueue
{
    struct multicastDelegate *events;
    struct eventQueue *next;  
};

struct delegateList
{
    void *target;
    void *sender;
    delegate del;
    struct delegateList *next;
};

struct multicastDelegate
{
    eventArgs *e;
    struct delegateList *subscribers;
    multicastDelegate *next;
};

eventQueue *rootQueue;
neuron *neurons[5];
int refCountF = 0;
float *NewFloat(float value)
{
    refCount++;
    refCountF++;
    float *f = (float*)malloc(sizeof(float));
    *f = value;
    return f;
}

float NoDelay(float signal)
{
	return 0.0f;
}

float InputDelay(float signal)
{
    float e = 2.718f;
    float a = 1.0f;
    float b = 0.75f;
    float c = 1.0f;
    float d = 2.0f;
    return c + (a / (pow(e, (b * signal) - d)));
}

float NoDecay(float potential, long deltaT)
{
    return potential;
}

float InstantDecay(float potential, long deltaT)
{
    if (deltaT > 0)
        return 0.0f;
    else
        return potential;
}

float ExponentialDecay(float potential, long deltaT)
{
    if (deltaT == 0 || potential == 0) return potential; // avoid the math to speed things up
    double dT = (double)deltaT;
    return (float)(potential / exp(0.25 * dT));
}

int refCountE = 0;
int _neurons = 5;

eventArgs* NewEventArgs(void* args, long t, dispose dispose)
{
    refCount++;
    refCountE++;
    eventArgs *e = (eventArgs*)malloc(sizeof(eventArgs));
    e->value = args;
    e->dispose = dispose;
    e->t = t;
    return e;
}

int refCountQ = 0;
eventQueue* NewQueue(eventQueue *parent)
{
    refCount++;
	refCountQ++;
    eventQueue* queue = (eventQueue*)malloc(sizeof(eventQueue));
    queue->events = 0;
    queue->next = 0;
    if (parent != 0) parent->next = queue;
    return queue;
}

int refCountM = 0;
multicastDelegate* NewMulticastDelegate(multicastDelegate *parent)
{
    refCount++;
	refCountM++;
    multicastDelegate* mcast = (multicastDelegate*)malloc(sizeof(multicastDelegate));
    mcast->e = 0;
    mcast->subscribers = 0;
    mcast->next = 0;
    if (parent != 0) parent->next = mcast;
    return mcast;
}

int refCountD = 0;
delegateList* NewDelegateList(delegateList *parent)
{
    refCount++;
	refCountD++;
    delegateList* del = (delegateList*)malloc(sizeof(delegateList));
    del->next = 0;
    del->target = 0;
    del->sender = 0;
    del->del = 0;
    if (parent != 0) parent->next = del;
    return del;
}

neuron* NewNeuron()
{
	refCount++;
	neuron *n = (neuron*)malloc(sizeof(neuron));
	n->threshold = 1.0f;
	n->potential = 0.0f;
	n->refracting = 0;
	n->Integration = 0;
	n->Activation = 0;
	n->Decay = 0;
	n->Delay = 0;
	n->axons = 0;
	n->dendrites = 0;
	n->isFiring = false;
	n->t = 0;
	return n;
}

int refCountS;
synapse* NewSynapse()
{
	refCount++;
	refCountS++;
	synapse *s = (synapse*)malloc(sizeof(synapse));
	s->input = 0;
	s->output = 0;
	s->weight = 1.0;
	return s;
}

int refCountSL = 0;
synapseList* NewSynapseList()
{
	refCount++;
	refCountSL++;
	synapseList *sl = (synapseList*)malloc(sizeof(synapseList));
	sl->next = 0;
	sl->value = 0;
	return sl;
}

multicastDelegate* NextMcastDelegate(eventQueue *queue)
{
    if (queue == 0) queue = NewQueue(0);
    multicastDelegate *mcast;
    if (queue->events == 0)
    {
		queue->events = NewMulticastDelegate(0);
		mcast = queue->events;
    }
    else
    {
		mcast = queue->events;
		while (mcast->next != 0)
		{
				mcast = mcast->next;
		}
		mcast = NewMulticastDelegate(mcast);
    }
    return mcast;
}

void FreeEventsArgs(eventArgs *e)
{
    if (e == 0) return;
    if (e->dispose != 0)
        e->dispose(e->value);
	_free((void **)&e);
    refCount--;
    refCountE--;
}

void FreeDelegateList(delegateList *dList)
{
    if (dList == 0) return;
    FreeDelegateList(dList->next);
	_free((void **)&dList);
    refCount--;
	refCountD--;
}

void FreeMulticast(multicastDelegate *mcast)
{
    if (mcast == 0) return;
	FreeMulticast(mcast->next);
    FreeDelegateList(mcast->subscribers);
	FreeEventsArgs(mcast->e);
	_free((void **)&mcast);
    refCount--;
	refCountM--;
}

void FreeQueue(eventQueue *queue)
{
    if (queue == 0) return;
    FreeMulticast(queue->events);
	_free((void **)&queue);
    refCount--;
	refCountQ--;
}

void FreeQueues(eventQueue *queue)
{
	if (queue == 0) return;
	FreeQueues(queue->next);
	FreeQueue(queue);
}

void FreeSynapse(synapse *s)
{
	if (s == 0) return;
	_free((void **)&s);
	refCount--;
	refCountS--;
}

synapseList* FreeSynapseList(synapseList *sl, bool isAxon)
{
	if (sl == 0) return sl;
	synapseList *rem = sl->next;
	if (isAxon) FreeSynapse(sl->value);
	_free((void **)&sl);
	refCount--;
	refCountSL--;
	return rem;
}

void FreeNeuron(neuron *n)
{
	if (n == 0) return;
	while (n->dendrites != 0)
		n->dendrites = FreeSynapseList(n->dendrites, false);
	while (n->axons != 0)
		n->axons = FreeSynapseList(n->axons, true);
	_free((void **)&n);
	refCount--;
}

void FreeNeurons()
{
	for (int i = 0; i < _neurons; i++)
	{
		FreeNeuron(neurons[i]);
	}
}

void FreeFloat(void* value)
{
	_free((void **)&value);
	refCount--;
	refCountF--;
}

void Integrate(int* target, int* sender, eventArgs *e)
{
	neuron *n = (struct neuron *)target;
	synapse *s = (synapse*)sender;
	float value = ((float*)e->value)[0] * s->weight;
	long t = e->t;
	long dT = t - n->t;
	n->t = t;
	n->potential = n->Decay(n->potential, dT);
	n->potential += value; // increment membrane potential
	if (n->potential >= n->threshold
			&& !n->isFiring)
	{
		// post a single delegate to the tail of the root queue to process the Activate
		// operation for the signaling neuron.
		// we need to defer the activation until after all the inputs have been processed
		// so that we can fully integrate all the input signals, prior to computing the 
		// activation.
		n->isFiring = true; // enforces single entry

		multicastDelegate *mcast = NextMcastDelegate(rootQueue);

		mcast->e = NewEventArgs(0, e->t, 0);

		delegateList *subscribers = NewDelegateList(0);
		subscribers->del = n->Activation; // calls the activation function
		subscribers->target = n;
		subscribers->sender = n;
		mcast->subscribers = subscribers;
	}
}

void L2Integrate(int* target, int* sender, eventArgs *e)
{
	Integrate(target, sender, e);
}

void EnqueueSpike(neuron *n, float value, long time)
{
	int delay = (int)n->Delay(n->potential / n->threshold);
	int t = delay;
	n->potential = 0.0;
	if (rootQueue == 0) rootQueue = NewQueue(0);
	eventQueue *queue = rootQueue;
	while (t > 1)
	{
		if (queue->next == 0)
		{
			// add a new queue
			queue = NewQueue(queue);
		}
		else queue = queue->next;
		t--;
	}
	// add a new event to the queue
	multicastDelegate *mcast = NextMcastDelegate(queue);
	mcast->e = NewEventArgs(NewFloat(value), time + delay, (dispose)FreeFloat);
	delegateList *subscribers = 0;
	delegateList *subs = 0;
	synapseList *axon = n->axons;

	while (axon != 0)
	{
		subs = NewDelegateList(subs);
		if (subscribers == 0) subscribers = subs;
		subs->del = axon->value->output->Integration;
		subs->target = axon->value->output;
		subs->sender = axon->value;
		axon = n->axons->next;
	}

	mcast->subscribers = subscribers;
}

void Activate(int* target, int* sender, eventArgs *e)
{
	long t = e->t;
	struct neuron *n = (struct neuron *)target;
	n->isFiring = false;
	if (n->potential >= n->threshold)
	{
		// spike!!!
		EnqueueSpike(n, 1.0f, e->t);
	}
}

void L2Activate(int* target, int* sender, eventArgs *e)
{
	Activate(target, sender, e);
}

void LoadNeurons()
{
	for (int i = 0; i < _neurons; i++)
	{
		neurons[i] = NewNeuron();
		neurons[i]->Integration = (delegate)Integrate;
		neurons[i]->Activation = (delegate)Activate;
		neurons[i]->Decay = (decay)InstantDecay;
		neurons[i]->Delay = (delay)InputDelay;
	}

	neurons[1]->Activation = (delegate)L2Activate;
	neurons[1]->Integration = (delegate)L2Integrate;
	neurons[1]->Decay = (decay)ExponentialDecay;
	neurons[1]->threshold = 3.0f;

	neurons[0]->axons = NewSynapseList();
	neurons[0]->axons->value = NewSynapse();
	neurons[0]->axons->value->input = neurons[0];
	neurons[0]->axons->value->output = neurons[1];
	neurons[0]->axons->value->weight = 1.0;

	neurons[1]->dendrites = NewSynapseList();
	neurons[1]->dendrites->value = neurons[0]->axons->value; // pointer to same synapse on n0 

	neurons[2]->axons = NewSynapseList();
	neurons[2]->axons->value = NewSynapse();
	neurons[2]->axons->value->input = neurons[0];
	neurons[2]->axons->value->output = neurons[1];
	neurons[2]->axons->value->weight = 1.0;

	neurons[1]->dendrites->next = NewSynapseList();
	neurons[1]->dendrites->next->value = neurons[2]->axons->value;
}

void ProcessDelegate(delegateList *delegateList, eventArgs *e)
{
	if (delegateList == 0) return;
	(*delegateList->del)(delegateList->target, delegateList->sender, e);
	ProcessDelegate(delegateList->next, e); // recurse the delegate list to invoke subsribers
}

void ProcessEvents(multicastDelegate *events, long t)
{
	if (events == 0) return;
	ProcessDelegate(events->subscribers, events->e);
	ProcessEvents(events->next, t);
}

void ProcessEventQueue(eventQueue *queue, long t)
{
	if (queue == 0) return;
	ProcessEvents(queue->events, t);
}

void ApplyInputs(long t)
{
	float v = 1.0;
	float *value = &v;
	eventArgs *e = NewEventArgs(NewFloat(1.0f), t, (dispose)FreeFloat);
	synapse *s = NewSynapse();
	s->weight = 1.0f;
	(*neurons[0]->Integration)((int*)neurons[0], (int*)s, e);
	FreeSynapse(s);
	FreeEventsArgs(e);
	//(*neurons[2]->Integration)((int*)neurons[2], (void*)value);
	//(*neurons[0]->Activation)((int*)neurons[0], (void*)0);
	//(*neurons[2]->Activation)((int*)neurons[2], (void*)0);
}

void Run(long t)
{
	ApplyInputs(t);
	ProcessEventQueue(rootQueue, t);
	eventQueue *rem = rootQueue->next;
	FreeQueue(rootQueue);
	rootQueue = rem;
}

int _tmain(int argc, _TCHAR* argv[])
{
	double t = 0;
	rootQueue = NewQueue(0);
	LoadNeurons();
	double ticks = 1000.0;
	double s = (double)clock() / ticks;
	while (t < 1000000)
	{
		Run((long)t++);
	}
	double e = (double)clock() / ticks;
	FreeQueues(rootQueue);
	FreeNeurons();
	double d = e - s;
	printf("Count: %f, Delta T: %f, Rate: %e, RefCount: %d", t, d, t / d, refCount);
	getchar();
	return 0;
}